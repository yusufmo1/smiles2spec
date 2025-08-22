"""Prediction API schemas."""
from pydantic import BaseModel, validator
from typing import List, Optional, Dict, Any

class PredictRequest(BaseModel):
    """Prediction request schema."""
    smiles: str
    
    @validator('smiles')
    def validate_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError('SMILES cannot be empty')
        # Basic SMILES validation
        cleaned = v.strip()
        if len(cleaned) > 1000:
            raise ValueError('SMILES too long (max 1000 characters)')
        return cleaned

class PeakSchema(BaseModel):
    """Peak schema."""
    mz: float
    intensity: float
    
    @validator('mz')
    def validate_mz(cls, v):
        if v < 0:
            raise ValueError('m/z must be non-negative')
        return v
    
    @validator('intensity')
    def validate_intensity(cls, v):
        if v < 0:
            raise ValueError('Intensity must be non-negative')
        return v

class SpectrumSchema(BaseModel):
    """Spectrum schema."""
    x: List[float]
    y: List[float]
    
    @validator('y')
    def validate_arrays_same_length(cls, v, values):
        if 'x' in values and len(v) != len(values['x']):
            raise ValueError('x and y arrays must have same length')
        return v

class PredictResponse(BaseModel):
    """Prediction response schema."""
    smiles: str
    chemical_name: str
    molecular_weight: float
    exact_mass: float
    spectrum: SpectrumSchema
    peaks: List[PeakSchema]
    structure_png: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None

class ExportMSPRequest(BaseModel):
    """MSP export request schema."""
    smiles: str
    
    @validator('smiles')
    def validate_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError('SMILES cannot be empty')
        return v.strip()

class ExportMSPBatchRequest(BaseModel):
    """Batch MSP export request schema."""
    smiles_list: List[str]
    
    @validator('smiles_list')
    def validate_smiles_list(cls, v):
        if not v:
            raise ValueError('SMILES list cannot be empty')
        if len(v) > 100:
            raise ValueError('Maximum 100 SMILES allowed per batch')
        
        # Validate each SMILES
        cleaned_list = []
        for smiles in v:
            if not smiles or not smiles.strip():
                continue  # Skip empty entries
            cleaned_list.append(smiles.strip())
        
        if not cleaned_list:
            raise ValueError('No valid SMILES found in list')
        
        return cleaned_list 