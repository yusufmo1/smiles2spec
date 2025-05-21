"""SMILES generation service implementation."""
import logging
from typing import List, Optional
from rdkit import Chem

from ..common.config import MODELS
from ..common.streaming import complete_response
from .prompts import SYSTEM_PROMPT, get_user_prompt
from .validator import extract_smiles_from_text, validate_smiles

logger = logging.getLogger(__name__)

def generate_smiles(
    description: str, 
    count: int = 1, 
    fallback: str = "C"
) -> List[str]:
    """Generate valid SMILES strings from a description."""
    try:
        # Create messages
        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": get_user_prompt(description, count)}
        ]
        
        # Get model response
        model = MODELS["chemistry"]
        logger.info(f"Requesting SMILES generation for: {description}")
        response_text = complete_response(messages, model)
        
        # Extract and validate SMILES
        smiles_candidates = extract_smiles_from_text(response_text)
        validated_results = validate_smiles(smiles_candidates)
        
        # Filter valid SMILES
        valid_smiles = [s for s, is_valid in validated_results if is_valid]
        
        if valid_smiles:
            # Limit to requested count
            return valid_smiles[:count]
        else:
            logger.warning(f"No valid SMILES generated for: {description}")
            return [fallback] * count
            
    except Exception as e:
        logger.error(f"SMILES generation error: {str(e)}")
        return [fallback] * count
