"""Spectrum domain models."""
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import numpy as np

@dataclass
class Peak:
    """Single spectrum peak."""
    mz: float
    intensity: float
    
    def to_dict(self) -> Dict[str, Any]:
        return {"mz": self.mz, "intensity": self.intensity}

@dataclass
class Spectrum:
    """Mass spectrum container."""
    x: np.ndarray  # m/z values
    y: np.ndarray  # intensity values
    peaks: List[Peak] = field(default_factory=list)
    
    def __post_init__(self):
        """Initialize peaks from x,y arrays if not provided."""
        if not self.peaks and len(self.x) == len(self.y):
            self.peaks = [Peak(mz=float(x), intensity=float(y)) 
                         for x, y in zip(self.x, self.y)]
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "x": self.x.tolist() if isinstance(self.x, np.ndarray) else self.x,
            "y": self.y.tolist() if isinstance(self.y, np.ndarray) else self.y,
            "peaks": [peak.to_dict() for peak in self.peaks]
        }

@dataclass
class PredictionResult:
    """Complete prediction result."""
    smiles: str
    chemical_name: str
    molecular_weight: float
    exact_mass: float
    spectrum: Spectrum
    structure_png: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to API response format."""
        result = {
            'smiles': self.smiles,
            'chemical_name': self.chemical_name,
            'molecular_weight': self.molecular_weight,
            'exact_mass': self.exact_mass,
            'spectrum': self.spectrum.to_dict(),
            'peaks': [peak.to_dict() for peak in self.spectrum.peaks]
        }
        
        if self.structure_png:
            result['structure_png'] = self.structure_png
        
        if self.metadata:
            result['metadata'] = self.metadata
        
        return result 