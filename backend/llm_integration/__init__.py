"""
LLM integration package for smiles2spec.
"""

# Import main functions for easy access
from .chat_service import generate_chat_response
from .smiles_generator import generate_random_smiles
from .image_analyzer import analyze_molecule_image
from ..model_downloader import initialize_models

__all__ = [
    'generate_chat_response',
    'generate_random_smiles',
    'analyze_molecule_image',
    'initialize_models'
] 