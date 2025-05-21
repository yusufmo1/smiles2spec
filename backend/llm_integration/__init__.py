"""
LLM integration package for the mass spectrometry prediction API.
"""

# Import main functions for easy access
from .smiles_generator.service import generate_smiles
from .spectra_chat.service import generate_chat_response

__all__ = [
    'generate_smiles',
    'generate_chat_response'
] 