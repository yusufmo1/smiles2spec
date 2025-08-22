"""LLM integration package."""
from .services.chat_service import ChatService
from .services.smiles_service import SMILESService

# Compatibility functions for existing imports
def generate_chat_response(messages, spectrum_data=None, stream=False):
    """Generate chat response (compatibility function)."""
    chat_service = ChatService()
    return chat_service.generate_response(messages, spectrum_data, stream)

def generate_smiles(count=1, description=""):
    """Generate SMILES strings (compatibility function)."""
    smiles_service = SMILESService()
    return smiles_service.generate_smiles(description, count)

__all__ = [
    'ChatService',
    'SMILESService', 
    'generate_chat_response',
    'generate_smiles'
] 