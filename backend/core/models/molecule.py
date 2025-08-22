"""Molecule domain models."""
from dataclasses import dataclass
from typing import Dict, List, Optional, Any
import numpy as np

@dataclass
class MolecularFeatures:
    """Molecular features container."""
    smiles: str
    descriptors: Optional[np.ndarray] = None
    descriptor_names: Optional[List[str]] = None
    fingerprints: Dict[str, np.ndarray] = None
    electronic_features: Optional[np.ndarray] = None
    atom_counts: Optional[Dict[str, int]] = None
    bond_counts: Optional[Dict[str, int]] = None
    
    def __post_init__(self):
        """Initialize default values."""
        if self.fingerprints is None:
            self.fingerprints = {}
        if self.atom_counts is None:
            self.atom_counts = {}
        if self.bond_counts is None:
            self.bond_counts = {}
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API responses."""
        result = {
            'smiles': self.smiles,
            'atom_counts': self.atom_counts,
            'bond_counts': self.bond_counts
        }
        
        if self.descriptors is not None:
            result['descriptors'] = self.descriptors.tolist()
            if self.descriptor_names:
                result['descriptor_names'] = self.descriptor_names
        
        if self.fingerprints:
            result['fingerprints'] = {
                name: fp.tolist() if isinstance(fp, np.ndarray) else fp 
                for name, fp in self.fingerprints.items()
            }
        
        if self.electronic_features is not None:
            result['electronic_features'] = self.electronic_features.tolist()
        
        return result

@dataclass
class MoleculeInfo:
    """Basic molecule information."""
    smiles: str
    chemical_name: str
    molecular_weight: float
    exact_mass: float
    structure_png: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API responses."""
        return {
            'smiles': self.smiles,
            'chemical_name': self.chemical_name,
            'molecular_weight': self.molecular_weight,
            'exact_mass': self.exact_mass,
            'structure_png': self.structure_png
        } 