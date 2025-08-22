"""Prediction service with clean interface."""
from typing import Optional, List
import os
import numpy as np

from ..config import settings
from ..core.models.molecule import MoleculeInfo
from ..core.models.spectrum import PredictionResult, Spectrum, Peak
from ..utils import logger, SMILESError, ModelError
from ..core.processors.feature_processor import FeatureProcessor
from ..core.processors.spectrum_processor import SpectralProcessor
from ..core.ml.model_handler import ModelHandler

class PredictionService:
    """
    High-level service for mass spectrum predictions.
    
    Orchestrates the complete prediction workflow from SMILES input
    to spectrum output. Manages component initialization, feature
    extraction, model prediction, and result formatting.
    
    This service provides a clean interface between the API layer
    and the core ML components, handling all complexity internally.
    
    Attributes:
        _model_handler: ML model handler for predictions
        _feature_processor: Molecular feature extractor
        _spectrum_processor: Spectrum data processor
        _initialized: Service initialization status
    """
    
    def __init__(self):
        self._model_handler = None
        self._feature_processor = None
        self._spectrum_processor = None
        self._initialized = False
    
    def initialize(self) -> bool:
        """
        Initialize all service components.
        
        Loads the ML model, preprocessor, and feature mapping from
        configured paths. Validates that required files exist before
        attempting initialization.
        
        Returns:
            bool: True if all components initialized successfully
            
        Note:
            Initialization is lazy - happens on first prediction if not
            explicitly called. Logs warnings for missing files.
        """
        try:
            # Check if model files exist
            model_path = settings.api.model_path
            preprocessor_path = settings.api.preprocessor_path
            feature_mapping_path = settings.api.feature_mapping_path
            
            if not os.path.exists(model_path):
                logger.warning(f"Model file not found: {model_path}")
                return False
            
            if not os.path.exists(preprocessor_path):
                logger.warning(f"Preprocessor file not found: {preprocessor_path}")
                return False
            
            # Initialize components
            self._feature_processor = FeatureProcessor()
            self._spectrum_processor = SpectralProcessor()
            self._model_handler = ModelHandler()
            
            # Load model and preprocessor
            success = self._model_handler.initialize(
                model_path, 
                preprocessor_path,
                feature_mapping_path if os.path.exists(feature_mapping_path) else None
            )
            
            if success:
                self._initialized = True
                logger.info("Prediction service initialized successfully")
            else:
                logger.error("Failed to initialize model handler")
                
            return success
            
        except Exception as e:
            logger.error(f"Error initializing prediction service: {str(e)}")
            return False
    
    def predict_spectrum(self, smiles: str) -> PredictionResult:
        """
        Predict mass spectrum for a single SMILES string.
        
        Complete prediction pipeline:
        1. Validate SMILES input
        2. Extract molecular features
        3. Generate structure information
        4. Predict spectrum intensities
        5. Format results with metadata
        
        Args:
            smiles: Valid SMILES molecular representation
            
        Returns:
            PredictionResult containing:
                - Predicted spectrum (m/z and intensities)
                - Molecular metadata (name, weight, formula)
                - Structure visualization (base64 PNG)
                
        Raises:
            SMILESError: Invalid SMILES input
            ModelError: Prediction service not initialized or prediction failed
            
        Note:
            Automatically initializes service on first call if needed
        """
        if not self._initialized:
            if not self.initialize():
                raise ModelError("Prediction service not initialized and initialization failed")
        
        try:
            # Validate and process SMILES
            if not smiles or not smiles.strip():
                raise SMILESError("Empty SMILES string provided")
            
            smiles = smiles.strip()
            
            # Extract molecular features
            mol_features = self._feature_processor.extract_features(smiles)
            mol_info = self._feature_processor.extract_molecule_info(smiles)
            
            # Make prediction
            if self._model_handler and self._model_handler.is_loaded():
                predicted_intensities = self._model_handler.predict(mol_features)
                
                # Convert to spectrum
                spectrum = self._spectrum_processor.create_spectrum_from_prediction(
                    predicted_intensities
                )
            else:
                # Create empty spectrum if model not available
                spectrum = self._spectrum_processor.create_empty_spectrum()
                logger.warning("Model not loaded, returning empty spectrum")
            
            # Create prediction result
            result = PredictionResult(
                smiles=mol_info.smiles,
                chemical_name=mol_info.chemical_name,
                molecular_weight=mol_info.molecular_weight,
                exact_mass=mol_info.exact_mass,
                spectrum=spectrum,
                structure_png=mol_info.structure_png
            )
            
            return result
            
        except SMILESError:
            raise
        except Exception as e:
            logger.error(f"Error predicting spectrum for SMILES {smiles}: {str(e)}")
            raise ModelError(f"Prediction failed: {str(e)}")
    
    def predict_batch(self, smiles_list: List[str]) -> List[PredictionResult]:
        """
        Predict spectra for multiple SMILES strings.
        
        Processes a batch of SMILES sequentially, collecting successful
        predictions and logging errors for failed ones. Useful for bulk
        processing and analysis.
        
        Args:
            smiles_list: List of SMILES strings to process
            
        Returns:
            List of PredictionResult objects for successful predictions
            
        Note:
            - Failed predictions are logged but don't stop batch processing
            - Empty list returned if all predictions fail
            - Consider implementing parallel processing for large batches
        """
        results = []
        
        for smiles in smiles_list:
            try:
                result = self.predict_spectrum(smiles)
                results.append(result)
            except Exception as e:
                logger.error(f"Error predicting spectrum for SMILES {smiles}: {str(e)}")
                # Could optionally include error results in batch
                continue
        
        return results
    
    def is_ready(self) -> bool:
        """
        Check if service is ready for predictions.
        
        Verifies that all required components are loaded and initialized.
        Used by health check endpoints to report service status.
        
        Returns:
            bool: True if service can make predictions
        """
        return self._initialized and self._model_handler and self._model_handler.is_loaded() 