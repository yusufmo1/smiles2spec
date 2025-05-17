"""
Main prediction service that integrates all components.
"""

import os
import numpy as np
from config import FEATURE_CONFIG, SPECTRAL_CONFIG, FEATURE_PROCESSING_CONFIG, API_CONFIG
from molecule_featurizer import compute_molecular_features
from spectrum_processor import SpectralProcessor
from model_handler import ModelHandler
from utils import logger, handle_error, convert_np_to_list, is_valid_smiles, smiles_to_png_base64, smiles_to_name

class PredictionService:
    def __init__(self):
        self.model_handler = ModelHandler()
        self.spectral_processor = SpectralProcessor(SPECTRAL_CONFIG)
        self.model_loaded = False
        self.initialize()
    
    def initialize(self):
        """Initialize the prediction service."""
        try:
            # Ensure model and processor files exist
            model_path = API_CONFIG.get('model_path', 'models/spectrum_predictor.pkl')
            preprocessor_path = API_CONFIG.get('preprocessor_path', 'models/feature_preprocessor.pkl')
            feature_mapping_path = API_CONFIG.get('feature_mapping_path', 'models/feature_mapping.json')
            
            if not os.path.exists(model_path):
                logger.warning(f"Model file not found: {model_path}")
                return False
            
            if not os.path.exists(preprocessor_path):
                logger.warning(f"Preprocessor file not found: {preprocessor_path}")
                return False
            
            # Initialize model handler
            self.model_loaded = self.model_handler.initialize(
                model_path, 
                preprocessor_path,
                feature_mapping_path if os.path.exists(feature_mapping_path) else None
            )
            
            return self.model_loaded
        except Exception as e:
            error_result = handle_error(e, "initializing prediction service")
            logger.error(f"Error initializing prediction service: {error_result['error']}")
            return False
    
    def predict_spectrum_from_smiles(self, smiles):
        """Predict mass spectrometry spectrum from SMILES."""
        try:
            # Validate SMILES
            if not is_valid_smiles(smiles):
                return {"error": f"Invalid SMILES string: {smiles}"}
            
            # Compute molecular features
            mol_features = compute_molecular_features(smiles, FEATURE_CONFIG)
            
            # Get chemical name
            chemical_name = smiles_to_name(smiles)
            
            # Make prediction if model is loaded
            if not self.model_loaded:
                # Return empty spectrum with molecular features
                molecular_weight, exact_mass = self._extract_molecular_weights(mol_features)
                return self._create_empty_response(smiles, molecular_weight, exact_mass, chemical_name)
            
            # Get prediction from model
            predicted_spectrum = self.model_handler.predict(mol_features)
            
            # Create response with binned spectrum
            mz_values = np.arange(0, len(predicted_spectrum) * SPECTRAL_CONFIG['bin_size'], SPECTRAL_CONFIG['bin_size'])
            peaks = self._convert_binned_to_peaks(mz_values, predicted_spectrum)
            
            # Extract molecular weights
            molecular_weight, exact_mass = self._extract_molecular_weights(mol_features)
            
            # Create response
            response = {
                "smiles": smiles,
                "chemical_name": chemical_name,
                "molecular_weight": molecular_weight,
                "exact_mass": exact_mass,
                "spectrum": {
                    "x": mz_values.tolist(),
                    "y": predicted_spectrum.tolist()
                },
                "peaks": peaks,
                "structure_png": smiles_to_png_base64(smiles)
            }
            
            return response
        except Exception as e:
            error_result = handle_error(e, f"predicting spectrum for SMILES {smiles}")
            return {"error": error_result['error'], "smiles": smiles}
    
    def _extract_molecular_weights(self, mol_features):
        """Extract molecular weight and exact mass from features."""
        molecular_weight = exact_mass = 0
        if 'descriptors' in mol_features and 'descriptor_names' in mol_features:
            try:
                molwt_idx = mol_features['descriptor_names'].index('MolWt')
                molecular_weight = float(mol_features['descriptors'][molwt_idx])
            except (ValueError, IndexError):
                pass
            try:
                exactwt_idx = mol_features['descriptor_names'].index('ExactMolWt')
                exact_mass = float(mol_features['descriptors'][exactwt_idx])
            except (ValueError, IndexError):
                pass
        return molecular_weight, exact_mass
    
    def _convert_binned_to_peaks(self, mz_values, intensities, threshold=0.01):
        """Convert binned spectrum to peak list, filtering by threshold."""
        peaks = []
        for mz, intensity in zip(mz_values, intensities):
            if intensity > threshold:  # Only include peaks above threshold
                peaks.append({"mz": float(mz), "intensity": float(intensity)})
        # Sort by intensity descending
        peaks.sort(key=lambda x: x["intensity"], reverse=True)
        return peaks
    
    def _create_empty_response(self, smiles, molecular_weight, exact_mass, chemical_name=""):
        """Create an empty response when model is not available."""
        # Create an empty binned spectrum
        binned_spectrum = np.zeros(int(SPECTRAL_CONFIG['max_mz'] / SPECTRAL_CONFIG['bin_size']) + 1)
        mz_values = np.arange(0, len(binned_spectrum) * SPECTRAL_CONFIG['bin_size'], SPECTRAL_CONFIG['bin_size'])
        
        return {
            "smiles": smiles,
            "chemical_name": chemical_name,
            "molecular_weight": molecular_weight,
            "exact_mass": exact_mass,
            "spectrum": {
                "x": mz_values.tolist(),
                "y": binned_spectrum.tolist()
            },
            "peaks": [],
            "warning": "No prediction model loaded. Returning empty spectrum with molecular information only.",
            "structure_png": smiles_to_png_base64(smiles)
        }