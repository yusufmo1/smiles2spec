"""API request handling utilities for LLM integration."""
import json
import requests
import logging
from typing import Dict, Any, Optional
from .config import OPENROUTER_BASE_URL, get_headers

logger = logging.getLogger(__name__)

def send_chat_request(
    messages: list, 
    model: str,
    temperature: float = 0.7,
    max_tokens: Optional[int] = None
) -> Dict[str, Any]:
    """
    Send a chat request to the OpenRouter API.
    
    Args:
        messages: List of message objects with role and content
        model: Model ID to use
        temperature: Sampling temperature (0-1)
        max_tokens: Maximum number of tokens to generate
        
    Returns:
        Full API response dictionary
    """
    url = f"{OPENROUTER_BASE_URL}/chat/completions"
    headers = get_headers()
    
    payload = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
    }
    
    if max_tokens:
        payload["max_tokens"] = max_tokens
    
    try:
        response = requests.post(url, headers=headers, json=payload)
        
        if response.status_code != 200:
            logger.error(f"API request failed: {response.status_code} - {response.text}")
            raise Exception(f"API request failed: {response.status_code} - {response.text}")
            
        return response.json()
        
    except Exception as e:
        logger.error(f"Error sending chat request: {str(e)}")
        raise

def extract_message_content(response: Dict[str, Any]) -> str:
    """Extract message content from API response."""
    if not response or "choices" not in response or not response["choices"]:
        raise ValueError("Invalid response format")
        
    return response["choices"][0]["message"]["content"]
