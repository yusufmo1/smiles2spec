"""Validation of generated SMILES strings."""
import re
from rdkit import Chem
from typing import List, Tuple


def validate_smiles(smiles_list: List[str]) -> List[Tuple[str, bool]]:
    """Validate a list of SMILES strings using RDKit."""
    results = []
    for smiles in smiles_list:
        cleaned = smiles.strip()
        if not cleaned:
            continue
        mol = Chem.MolFromSmiles(cleaned)
        is_valid = mol is not None
        results.append((cleaned, is_valid))
    return results


def extract_smiles_from_text(text: str) -> List[str]:
    """Extract SMILES strings from the model's response text."""
    lines = [line.strip() for line in text.split("\n")]
    if '```' in text:
        code_blocks = re.findall(r'```(?:smiles)?\n(.*?)```', text, re.DOTALL)
        if code_blocks:
            lines = []
            for block in code_blocks:
                lines.extend([line.strip() for line in block.split('\n')])
    cleaned_lines = []
    for line in lines:
        if not line:
            continue
        cleaned = re.sub(r'^\d+[\.\)]\s*', '', line)
        cleaned = cleaned.replace('`', '')
        cleaned_lines.append(cleaned)
    return cleaned_lines
