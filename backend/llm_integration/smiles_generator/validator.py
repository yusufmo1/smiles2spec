"""Validation of generated SMILES strings."""
import re
import logging
from rdkit import Chem
from typing import List, Tuple

logger = logging.getLogger(__name__)

def validate_smiles(smiles_list: List[str]) -> List[Tuple[str, bool]]:
    """Validate a list of SMILES strings using RDKit."""
    results = []
    
    for smiles in smiles_list:
        # Basic cleanup
        cleaned = smiles.strip()
        
        # Skip empty lines
        if not cleaned:
            continue
            
        # Validate with RDKit
        mol = Chem.MolFromSmiles(cleaned)
        is_valid = mol is not None
        
        if is_valid:
            logger.info(f"Valid SMILES: {cleaned}")
        else:
            logger.warning(f"Invalid SMILES: {cleaned}")
        
        results.append((cleaned, is_valid))
        
    return results

def extract_smiles_from_text(text: str) -> List[str]:
    """Extract SMILES strings from the model's response text."""
    logger.debug(f"Extracting SMILES from text: {text[:100]}...")
    
    # List to store potential SMILES
    candidates = []
    
    # Split by lines and clean up
    lines = [line.strip() for line in text.split('\n')]
    
    # Try to handle code blocks specifically
    if '```' in text:
        # Extract content from code blocks
        code_blocks = re.findall(r'```(?:smiles)?\n(.*?)```', text, re.DOTALL)
        if code_blocks:
            logger.debug(f"Found code blocks: {len(code_blocks)}")
            for block in code_blocks:
                candidates.extend([line.strip() for line in block.split('\n') if line.strip()])
    
    # If no candidates from code blocks, process line by line
    if not candidates:
        for line in lines:
            # Skip empty lines
            if not line:
                continue
                
            # Remove numbering like "1. " or "1) "
            cleaned = re.sub(r'^\d+[\.\)]\s*', '', line)
            # Remove backticks if any
            cleaned = cleaned.replace('`', '')
            
            # Skip if it's just explanatory text (heuristic check)
            if len(cleaned) > 5 and ' ' in cleaned and len(cleaned.split()) > 3:
                continue
                
            candidates.append(cleaned)
    
    # Additional cleaning: some models output things like "SMILES: CCO"
    final_candidates = []
    for candidate in candidates:
        # Extract SMILES after a label like "SMILES:" or "Structure:"
        match = re.search(r'(?:SMILES|Structure|Molecule|Compound):\s*([^\s]+)', candidate)
        if match:
            final_candidates.append(match.group(1))
        else:
            final_candidates.append(candidate)
    
    logger.debug(f"Extracted candidates: {final_candidates}")
    return final_candidates
