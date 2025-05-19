"""
SMILES generation service using LLM approach.
"""

from rdkit import Chem
import os
import json
import requests
import re
from ..utils import logger
from dotenv import load_dotenv
from .llm_config import (
    get_model_for_task,
    get_system_prompt,
    get_openrouter_headers,
    OPENROUTER_BASE_URL
)

# Load environment variables from .env file
dotenv_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(dotenv_path)

def generate_random_smiles(count=1, description="general organic molecule"):
    """
    Generate SMILES strings using LLM.
    
    Args:
        count: Number of SMILES strings to generate (default: 1)
        description: Text description of the desired molecule type (default: "general organic molecule")
        
    Returns:
        List of generated SMILES strings
    """
    try:
        # Generate SMILES using the LLM
        llm_smiles = generate_smiles_from_description(description, count)
        if llm_smiles:
            return llm_smiles
                
        # If LLM generation fails, return a simple molecule
        logger.warning(f"LLM SMILES generation failed for: {description}")
        return ["C"] * count  # Methane as ultimate fallback
    except Exception as e:
        logger.error(f"SMILES generation error: {str(e)}")
        # Return a simple molecule if generation fails
        return ["C"] * count  # Methane as ultimate fallback

def generate_smiles_from_description(description, count=1):
    """
    Generate SMILES strings based on text description using LLM.
    
    Args:
        description: Text description of desired molecules
        count: Number of SMILES to generate
        
    Returns:
        List of SMILES strings or None if generation failed
    """
    try:
        # Check if API key is available
        api_key = os.environ.get('OPENROUTER_API_KEY')
        if not api_key:
            logger.error("OpenRouter API key missing")
            return None
        
        # Get model and headers
        model = get_model_for_task("chemistry")
        headers = get_openrouter_headers()
        system_prompt = get_system_prompt("smiles_generator")
        
        # Request data
        data = {
            "model": model,
            "messages": [
                {
                    "role": "system",
                    "content": system_prompt
                },
                {
                    "role": "user",
                    "content": f"Generate {count} chemically valid SMILES strings for: {description}"
                }
            ]
        }
        
        # Make API request
        logger.info(f"Requesting SMILES from LLM for: {description}")
        response = requests.post(
            f"{OPENROUTER_BASE_URL}/chat/completions",
            headers=headers,
            data=json.dumps(data)
        )
        
        if response.status_code != 200:
            logger.error(f"API request failed: {response.text}")
            return None
            
        result = response.json()
        
        # Check for errors
        if "error" in result:
            logger.error(f"API returned error: {result['error'].get('message', 'Unknown error')}")
            if "metadata" in result.get("error", {}):
                logger.error(f"Error details: {result['error']['metadata']}")
            return None
            
        # Extract content
        if "choices" in result and len(result["choices"]) > 0:
            content = result["choices"][0]["message"]["content"]
            return process_smiles_response(content, count)
        else:
            logger.error("No choices found in API response")
            return None
        
    except Exception as e:
        logger.error(f"LLM SMILES generation error: {str(e)}")
        return None

def process_smiles_response(content, count):
    """Process LLM response and extract valid SMILES"""
    # Parse SMILES strings from the response
    smiles_list = [line.strip() for line in content.split('\n') if line.strip()]
    
    # Clean up SMILES (remove markdown code blocks, numbering, etc.)
    cleaned_smiles = []
    for s in smiles_list:
        # Remove code formatting
        s = s.replace('`', '').strip()
        
        # Remove numbering patterns like "1. " or "1) " at the beginning of lines
        s = re.sub(r'^\d+[\.\)]\s*', '', s)
        
        if s:
            cleaned_smiles.append(s)
    
    # Validate each SMILES
    valid_smiles = []
    for s in cleaned_smiles:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            valid_smiles.append(s)
            logger.info(f"Valid SMILES generated: {s}")
        else:
            logger.warning(f"Invalid SMILES generated: {s}")
    
    # Return valid SMILES up to the requested count
    return valid_smiles[:count] if valid_smiles else None 