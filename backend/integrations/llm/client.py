"""LLM client abstraction."""
from typing import List, Dict, Any, Generator, Optional
import requests
import json
from ...config import settings
from ...utils import logger, ExternalServiceError

class LLMClient:
    """
    Client for interacting with Large Language Models via OpenRouter API.
    
    Provides a unified interface for chat completions with various LLMs
    including GPT-4, google/gemma-3-27b-it:free, and others available through OpenRouter.
    Supports both streaming and non-streaming responses.
    
    OpenRouter acts as a gateway to multiple LLM providers, simplifying
    API management and providing fallback options.
    
    Attributes:
        config: LLM configuration from settings
        base_url: OpenRouter API endpoint
    """
    
    def __init__(self):
        self.config = settings.llm
        self.base_url = "https://openrouter.ai/api/v1"
    
    def _get_headers(self) -> Dict[str, str]:
        """
        Get required headers for OpenRouter API requests.
        
        Returns:
            Dictionary with authorization and metadata headers
            
        Raises:
            ExternalServiceError: If API key not configured
            
        Headers include:
            - Authorization: Bearer token
            - HTTP-Referer: Application URL for tracking
            - X-Title: Application name for OpenRouter dashboard
        """
        if not self.config.openrouter_api_key:
            raise ExternalServiceError("OPENROUTER_API_KEY not set", "openrouter")
        
        return {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.config.openrouter_api_key}",
            "HTTP-Referer": self.config.site_url,
            "X-Title": self.config.site_name,
        }
    
    def chat_completion(
        self, 
        messages: List[Dict[str, str]], 
        model: str,
        stream: bool = False,
        temperature: float = 0.7
    ) -> Any:
        """
        Send chat completion request to LLM.
        
        Args:
            messages: List of message dicts with 'role' and 'content'
            model: Model identifier (e.g., 'openai/gpt-4')
            stream: Whether to stream response tokens
            temperature: Sampling temperature (0-2, higher = more random)
            
        Returns:
            str: Complete response (non-streaming)
            Generator[str]: Token generator (streaming)
            
        Raises:
            ExternalServiceError: API request failed
            
        Example:
            messages = [
                {"role": "system", "content": "You are a chemist."},
                {"role": "user", "content": "What is benzene?"}
            ]
            response = client.chat_completion(messages, "openai/gpt-4")
        """
        url = f"{self.base_url}/chat/completions"
        headers = self._get_headers()
        
        payload = {
            "model": model,
            "messages": messages,
            "temperature": temperature,
            "stream": stream
        }
        
        try:
            if stream:
                return self._stream_request(url, headers, payload)
            else:
                response = requests.post(url, headers=headers, json=payload, timeout=30)
                if response.status_code != 200:
                    raise ExternalServiceError(f"API request failed: {response.status_code}", "openrouter")
                return response.json()["choices"][0]["message"]["content"]
        except Exception as e:
            logger.error(f"LLM API error: {str(e)}")
            raise ExternalServiceError(f"LLM request failed: {str(e)}", "openrouter")
    
    def _stream_request(self, url: str, headers: Dict, payload: Dict) -> Generator[str, None, None]:
        """
        Handle server-sent event (SSE) streaming responses.
        
        Processes streaming responses from OpenRouter, yielding content
        tokens as they arrive. Handles buffering and parsing of SSE format.
        
        Args:
            url: API endpoint URL
            headers: Request headers
            payload: Request payload
            
        Yields:
            str: Individual content tokens from the response
            
        Note:
            - Handles incomplete chunks and buffering
            - Ignores empty deltas and parse errors
            - Yields error messages on failure
        """
        try:
            with requests.post(url, headers=headers, json=payload, stream=True, timeout=30) as r:
                if r.status_code != 200:
                    raise ExternalServiceError(f"Stream request failed: {r.status_code}", "openrouter")
                
                buffer = ""
                for chunk in r.iter_content(chunk_size=1024, decode_unicode=True):
                    if not chunk:
                        continue
                    
                    buffer += chunk
                    while '\n' in buffer:
                        line, buffer = buffer.split('\n', 1)
                        line = line.strip()
                        
                        if line.startswith('data: '):
                            data = line[6:]
                            if data == '[DONE]':
                                return
                            
                            try:
                                data_obj = json.loads(data)
                                content = data_obj["choices"][0]["delta"].get("content", "")
                                if content:
                                    yield content
                            except json.JSONDecodeError:
                                continue
        except Exception as e:
            logger.error(f"Streaming error: {str(e)}")
            yield f"Error: {str(e)}" 