"""Chemistry utility functions."""
import base64
import re
from typing import Optional
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolDescriptors

from .logging import logger

def is_valid_smiles(smi: str) -> bool:
    """Check if SMILES string is valid."""
    if not isinstance(smi, str):
        return False
    try:
        mol = Chem.MolFromSmiles(smi)
        return mol is not None
    except Exception as e:
        logger.error(f"Error validating SMILES {smi}: {str(e)}")
        return False

def smiles_to_png_base64(smiles: str, w: int = 960, h: int = 720) -> str:
    """Convert SMILES string to PNG representation of the molecule encoded as base64.
    
    Args:
        smiles: SMILES string of the molecule
        w: Width of the PNG image
        h: Height of the PNG image
        
    Returns:
        Base64-encoded PNG image of the molecule
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        
        drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
        opts = drawer.drawOptions()
        opts.clearBackground = False  # transparent alpha channel
        
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        
        png_bytes = drawer.GetDrawingText()
        return base64.b64encode(png_bytes).decode()
    except Exception as e:
        logger.error(f"Error generating PNG for SMILES {smiles}: {str(e)}")
        return ""

def smiles_to_name(smiles: str) -> str:
    """Return IUPAC name if RDKit can, else empirical formula."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        
        # Try IUPAC (works if RDKit was compiled with NCI resolver)
        try:
            name = Chem.MolToIUPACName(mol)
            if name:
                return name
        except Exception:
            pass
        
        # Fallback to formula
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        logger.error(f"Error getting name for SMILES {smiles}: {str(e)}")
        return ""

def safe_filename(smiles: str) -> str:
    """Create a Windows / *nix friendly filename from a SMILES string."""
    return re.sub(r'[^A-Za-z0-9\-_\.]+', '_', smiles)[:120] or "spectrum"

def smiles_to_formula(smiles: str) -> str:
    """Convert SMILES to molecular formula."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return rdMolDescriptors.CalcMolFormula(mol)
        return ""
    except Exception as e:
        logger.error(f"Error calculating formula for SMILES {smiles}: {str(e)}")
        return "" 