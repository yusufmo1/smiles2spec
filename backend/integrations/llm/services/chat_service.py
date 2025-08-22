"""Chat service implementation."""
from typing import List, Dict, Any, Optional, Generator, Union
from ..client import LLMClient
from ....core.models.spectrum import PredictionResult
from ....utils import logger

class ChatService:
    """
    Chat service for AI-powered conversations about mass spectrometry.
    
    Provides context-aware chat functionality where the AI assistant
    can reference and discuss specific mass spectra and molecular
    structures. Uses google/gemma-3-27b-it:free for high-quality responses.
    
    Features:
    - Context injection for spectrum-specific discussions
    - Streaming and non-streaming response modes
    - Error handling with graceful fallbacks
    - Scientific accuracy focus for chemistry topics
    
    Attributes:
        client: LLMClient instance for API communication
    """
    
    def __init__(self):
        self.client = LLMClient()
    
    def generate_response(
        self,
        messages: List[Dict[str, str]],
        spectrum_context: Optional[Union[PredictionResult, Dict[str, Any]]] = None,
        stream: bool = False
    ) -> Any:
        """
        Generate AI response with optional mass spectrum context.
        
        Enhances the conversation with specific molecular and spectral
        information when available, allowing the AI to provide targeted
        insights about the current analysis.
        
        Args:
            messages: Chat history with role/content pairs
            spectrum_context: Optional spectrum data (PredictionResult or dict)
            stream: Whether to stream the response
            
        Returns:
            str: Complete response (non-streaming)
            Generator[str]: Token stream (streaming)
            
        Note:
            - Automatically prepends system prompt with context
            - Handles errors gracefully with error messages
            - Uses google/gemma-3-27b-it:free for best chemistry understanding
        """
        try:
            # Add system prompt with context
            system_prompt = self._build_system_prompt(spectrum_context)
            full_messages = [{"role": "system", "content": system_prompt}] + messages
            
            return self.client.chat_completion(
                messages=full_messages,
                model="google/gemma-3-27b-it:free",
                stream=stream
            )
        except Exception as e:
            logger.error(f"Chat service error: {str(e)}")
            if stream:
                def error_gen():
                    yield f"Error: {str(e)}"
                return error_gen()
            else:
                return f"Error: {str(e)}"
    
    def _build_system_prompt(self, spectrum_context: Optional[Union[PredictionResult, Dict[str, Any]]]) -> str:
        """
        Build system prompt with optional spectrum context.
        
        Creates a comprehensive system prompt that establishes the AI's
        role as a mass spectrometry expert and includes specific molecular
        context when available.
        
        Args:
            spectrum_context: Optional spectrum/molecule data
            
        Returns:
            str: Complete system prompt with or without context
            
        Context includes:
            - SMILES notation
            - Chemical name and formula
            - Molecular weight and exact mass
            - Number of spectral peaks
            
        Note:
            Handles both PredictionResult objects and dictionary formats
            for flexibility with different data sources.
        """
        base_prompt = """You are Spectra, a helpful AI assistant specializing in mass spectrometry and molecular chemistry.

You can help with:
- Interpreting mass spectra and fragmentation patterns
- Explaining SMILES notation and chemical structures
- Analyzing molecular properties and relationships
- Providing insights into chemical analysis workflows

Please provide clear, scientifically accurate responses."""
        
        if spectrum_context:
            # Handle both PredictionResult objects and dicts
            if hasattr(spectrum_context, 'smiles'):
                # PredictionResult object
                context = f"""

Currently analyzing:
- SMILES: {spectrum_context.smiles}
- Chemical Name: {spectrum_context.chemical_name}
- Molecular Weight: {spectrum_context.molecular_weight:.2f} g/mol
- Exact Mass: {spectrum_context.exact_mass:.4f} amu
- Number of Peaks: {len(spectrum_context.spectrum.peaks)}
"""
            else:
                # Dict format (from .to_dict())
                context = f"""

Currently analyzing:
- SMILES: {spectrum_context.get('smiles', 'Unknown')}
- Chemical Name: {spectrum_context.get('chemical_name', 'Unknown')}
- Molecular Weight: {spectrum_context.get('molecular_weight', 0):.2f} g/mol
- Exact Mass: {spectrum_context.get('exact_mass', 0):.4f} amu
- Number of Peaks: {len(spectrum_context.get('peaks', []))}
"""
            return base_prompt + context + "\n\nUse this context to provide specific insights about this molecule when relevant."
        
        return base_prompt 