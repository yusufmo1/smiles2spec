"""Validation of generated SMILES strings."""
import re
from rdkit import Chem
from typing import List, Tuple


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

        results.append((cleaned, is_valid))

    return results


def extract_smiles_from_text(text: str) -> List[str]:
    """Extract SMILES strings from the model's response text."""
    # Split by lines and clean up
    lines = [line.strip() for line in text.split('\n')]

    # Remove code blocks markers if present
    if '```' in text:
        # Extract content from code blocks
        code_blocks = re.findall(r'```(?:smiles)?\n(.*?)```', text, re.DOTALL)
        if code_blocks:
            lines = []
            for block in code_blocks:
                lines.extend([line.strip() for line in block.split('\n')])

    # Remove numbering patterns like "1. " or "1) "
    cleaned_lines = []
    for line in lines:
        # Skip empty lines
        if not line:
            continue

        # Remove numbering
        cleaned = re.sub(r'^\d+[\.\)]\s*', '', line)
        # Remove backticks if any
        cleaned = cleaned.replace('`', '')

        cleaned_lines.append(cleaned)

    return cleaned_lines
