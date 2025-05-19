"""
SMILES generation service using both deterministic and LLM approaches.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import random
import os
import requests
from ..utils import logger

def generate_random_smiles(count=1, description=None):
    """
    Generate random SMILES strings.
    
    Args:
        count: Number of SMILES strings to generate (default: 1)
        description: Text description of the desired molecule type (default: None)
        
    Returns:
        List of generated SMILES strings
    """
    try:
        # If description is provided, try to use LLM approach
        if description and description.strip():
            llm_smiles = generate_smiles_from_description(description, count)
            if llm_smiles:
                return llm_smiles
                
        # Fall back to deterministic generation
        return generate_deterministic_smiles(count)
    except Exception as e:
        logger.error(f"SMILES generation error: {str(e)}")
        # Return a simple molecule if generation fails
        return ["C"] * count  # Methane as ultimate fallback

def generate_deterministic_smiles(count=1):
    """
    Generate random SMILES using deterministic approach with RDKit.
    
    Args:
        count: Number of SMILES strings to generate
        
    Returns:
        List of generated SMILES strings
    """
    # Define common molecular fragments for medicinal chemistry
    fragments = [
        "c1ccccc1", "C1CCCCC1", "c1ccncc1", "C1CCNCC1", 
        "CC(=O)O", "CCO", "CN", "CF", "CCl", "CBr", "NC=O", 
        "C(=O)O", "c1cc(F)ccc1", "c1cc(Cl)ccc1", "CC#N", "C=C",
        "CCOCC", "CN(C)C", "CSC", "CC(=O)N", "NC(=O)O", "COC",
        "c1ccc(O)cc1", "c1ccc(N)cc1", "c1nc[nH]c1", "c1cnoc1"
    ]
    
    result = []
    for _ in range(count):
        # Generate a random SMILES by combining fragments
        n_fragments = random.randint(1, 3)
        selected = random.sample(fragments, n_fragments)
        
        # Create molecule
        mol = None
        for i, frag in enumerate(selected):
            frag_mol = Chem.MolFromSmiles(frag)
            if i == 0:
                mol = frag_mol
            else:
                # Simple merge (this is a simplified approach)
                combo = Chem.CombineMols(mol, frag_mol)
                mol = combo
        
        # If random selection failed, return a default SMILES
        if mol is None:
            smiles = "CCO"  # Ethanol as fallback
        else:
            smiles = Chem.MolToSmiles(mol)
            
        result.append(smiles)
    
    return result

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
        # This is where you would integrate with an LLM API
        # For example with OpenAI API:
        
        # api_key = os.environ.get('OPENAI_API_KEY')
        # if not api_key:
        #    return None
        #
        # response = requests.post(
        #     "https://api.openai.com/v1/chat/completions",
        #     headers={
        #         "Content-Type": "application/json",
        #         "Authorization": f"Bearer {api_key}"
        #     },
        #     json={
        #         "model": "gpt-4",
        #         "messages": [
        #             {"role": "system", "content": 
        #              "You are a chemical structure generator. Generate valid SMILES strings based on descriptions. "
        #              "Only output the SMILES strings, one per line, nothing else."},
        #             {"role": "user", "content": 
        #              f"Generate {count} chemically valid SMILES strings for: {description}"}
        #         ],
        #         "temperature": 0.7
        #     }
        # )
        #
        # if response.status_code == 200:
        #     result = response.json()
        #     content = result["choices"][0]["message"]["content"]
        #     # Parse SMILES strings from the response
        #     smiles_list = [line.strip() for line in content.split('\n') if line.strip()]
        #     
        #     # Validate each SMILES
        #     valid_smiles = []
        #     for s in smiles_list:
        #         mol = Chem.MolFromSmiles(s)
        #         if mol is not None:
        #             valid_smiles.append(s)
        #     
        #     # Return valid SMILES up to the requested count
        #     return valid_smiles[:count] if valid_smiles else None
        
        # For now, return None to fall back to deterministic generation
        return None
        
    except Exception as e:
        logger.error(f"LLM SMILES generation error: {str(e)}")
        return None 