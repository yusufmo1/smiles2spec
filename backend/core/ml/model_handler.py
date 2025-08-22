"""Model loading and prediction functionality."""
import os
import numpy as np
import json
from typing import Dict, Optional
from joblib import load

from ...config import settings
from ...core.models.molecule import MolecularFeatures
from ...utils import logger, ModelError, ensure_numpy_array

class ModelHandler:
    """
    Handler for machine learning model operations.
    
    Manages the loading, initialization, and prediction workflow for the
    mass spectrometry prediction model. Handles feature preprocessing,
    model inference, and feature mapping for consistent predictions.
    
    Attributes:
        model: Loaded scikit-learn model (Random Forest)
        preprocessor: Feature preprocessor for scaling and normalization
        feature_mapping: Optional mapping for consistent feature ordering
        _loaded: Boolean flag indicating if all components are loaded
    """
    
    def __init__(self):
        self.model = None
        self.preprocessor = None
        self.feature_mapping = None
        self._loaded = False
    
    def initialize(self, model_path: str, preprocessor_path: str, mapping_path: Optional[str] = None) -> bool:
        """
        Initialize the model handler with all required components.
        
        Loads the trained model, preprocessor, and optional feature mapping.
        All components must load successfully for initialization to complete.
        
        Args:
            model_path: Path to the saved scikit-learn model (.joblib)
            preprocessor_path: Path to the saved preprocessor (.joblib)
            mapping_path: Optional path to feature mapping JSON file
            
        Returns:
            bool: True if all components loaded successfully, False otherwise
            
        Note:
            The model and preprocessor are required; feature mapping is optional
            but recommended for production deployments.
        """
        try:
            model_loaded = self._load_model(model_path)
            preprocessor_loaded = self._load_preprocessor(preprocessor_path)
            
            mapping_loaded = True
            if mapping_path:
                mapping_loaded = self._load_feature_mapping(mapping_path)
            
            self._loaded = model_loaded and preprocessor_loaded and mapping_loaded
            
            if self._loaded:
                logger.info("Model handler initialized successfully")
            else:
                logger.error("Failed to initialize model handler")
                
            return self._loaded
            
        except Exception as e:
            logger.error(f"Error initializing model handler: {str(e)}")
            return False
    
    def predict(self, mol_features: MolecularFeatures) -> np.ndarray:
        """
        Make a mass spectrum prediction using the loaded model.
        
        Takes molecular features and predicts the mass spectrum intensities.
        Features are preprocessed and aligned according to the model's training.
        
        Args:
            mol_features: MolecularFeatures object containing:
                - descriptors: RDKit molecular descriptors
                - fingerprints: Various molecular fingerprints
                - electronic_features: Electronic and quantum properties
                
        Returns:
            np.ndarray: Predicted intensities for mass spectrum (1D array)
            
        Raises:
            ModelError: If model not loaded or prediction fails
        """
        if not self._loaded or self.model is None:
            raise ModelError("Model not loaded. Call initialize() first.")
        
        try:
            # Prepare features
            X = self._prepare_features(mol_features)
            
            # Make prediction
            prediction = self.model.predict(X)
            
            # Return as 1D array
            return prediction[0]
            
        except Exception as e:
            logger.error(f"Error making prediction: {str(e)}")
            raise ModelError(f"Prediction failed: {str(e)}")
    
    def is_loaded(self) -> bool:
        """
        Check if model is loaded and ready for predictions.
        
        Returns:
            bool: True if model and all required components are loaded
        """
        return self._loaded and self.model is not None
    
    def _load_model(self, model_path: str) -> bool:
        """
        Load a trained scikit-learn model from file.
        
        Args:
            model_path: Path to the model file (.joblib format)
            
        Returns:
            bool: True if model loaded successfully, False otherwise
            
        Note:
            Uses joblib for efficient loading of scikit-learn models
        """
        try:
            if not os.path.exists(model_path):
                logger.error(f"Model file not found: {model_path}")
                return False
            
            self.model = load(model_path)
            logger.info(f"Loaded model from {model_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error loading model from {model_path}: {str(e)}")
            return False
    
    def _load_preprocessor(self, preprocessor_path: str) -> bool:
        """
        Load a feature preprocessor from file.
        
        The preprocessor handles feature scaling, normalization, and
        missing value imputation to match the model's training data.
        
        Args:
            preprocessor_path: Path to the preprocessor file (.joblib)
            
        Returns:
            bool: True if preprocessor loaded successfully, False otherwise
        """
        try:
            if not os.path.exists(preprocessor_path):
                logger.error(f"Preprocessor file not found: {preprocessor_path}")
                return False
            
            # Import here to avoid circular imports
            from ..processors.feature_preprocessor import FeaturePreprocessor
            self.preprocessor = FeaturePreprocessor.load(preprocessor_path)
            logger.info(f"Loaded preprocessor from {preprocessor_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error loading preprocessor from {preprocessor_path}: {str(e)}")
            return False
    
    def _load_feature_mapping(self, mapping_path: str) -> bool:
        """
        Load feature mapping configuration from JSON file.
        
        The mapping ensures features are arranged in the exact order
        expected by the trained model, critical for accurate predictions.
        
        Args:
            mapping_path: Path to the feature mapping JSON file
            
        Returns:
            bool: True if mapping loaded successfully, False otherwise
            
        JSON Structure:
            - features: Dict mapping feature names to indices
            - metadata: Information about feature dimensions
        """
        try:
            if not os.path.exists(mapping_path):
                logger.error(f"Feature mapping file not found: {mapping_path}")
                return False
            
            with open(mapping_path, 'r') as f:
                self.feature_mapping = json.load(f)
            logger.info(f"Loaded feature mapping from {mapping_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error loading feature mapping from {mapping_path}: {str(e)}")
            return False
    
    def _prepare_features(self, mol_features: MolecularFeatures) -> np.ndarray:
        """
        Extract and preprocess features from molecular features.
        
        Intelligently selects between mapped and unmapped feature preparation
        based on available configuration. Ensures features are properly
        formatted for model input.
        
        Args:
            mol_features: MolecularFeatures object with computed features
            
        Returns:
            np.ndarray: Preprocessed feature array ready for prediction
            
        Raises:
            ModelError: If feature preparation fails
        """
        try:
            # If we have feature mapping, use structured approach
            if self.feature_mapping and (
                'filtered_names' in self.feature_mapping or
                'feature_order' in self.feature_mapping):
                return self._prepare_features_with_mapping(mol_features)
            else:
                # Fallback to original method if no mapping available
                return self._prepare_features_without_mapping(mol_features)
                
        except Exception as e:
            logger.error(f"Error preparing features: {str(e)}")
            raise ModelError(f"Feature preparation failed: {str(e)}")
    
    def _prepare_features_without_mapping(self, mol_features: MolecularFeatures) -> np.ndarray:
        """
        Prepare features using predetermined concatenation order.
        
        Used when no feature mapping is available. Concatenates features
        in a predetermined order that must match the model's training.
        
        Args:
            mol_features: MolecularFeatures object
            
        Returns:
            np.ndarray: Feature array with shape (1, n_features)
            
        Feature Order:
            1. RDKit descriptors
            2. Morgan fingerprint
            3. Morgan feature fingerprint
            4. Atom pairs
            5. RDKit fingerprint
            6. Avalon fingerprint
            7. Pattern fingerprint
            8. Layered fingerprint
            9. Electronic features
        """
        features = []
        
        # Add descriptors
        if mol_features.descriptors is not None:
            features.append(ensure_numpy_array(mol_features.descriptors))
        
        # Add fingerprints in expected order
        fingerprint_keys = [
            'morgan_fingerprint', 'morgan_feature_fp', 'atom_pairs', 
            'rdkit_fingerprint', 'avalon_fingerprint', 'pattern_fingerprint', 
            'layered_fingerprint'
        ]
        
        for key in fingerprint_keys:
            if key in mol_features.fingerprints and mol_features.fingerprints[key] is not None:
                features.append(ensure_numpy_array(mol_features.fingerprints[key]))
        
        # Add electronic features
        if mol_features.electronic_features is not None:
            features.append(ensure_numpy_array(mol_features.electronic_features))
        
        if not features:
            raise ModelError("No valid features found in molecular data")
        
        # Concatenate features
        X = np.concatenate(features).reshape(1, -1)
        
        # Preprocess if preprocessor is available
        if self.preprocessor:
            X = self.preprocessor.transform(X)
        
        return X
    
    def _prepare_features_with_mapping(self, mol_features: MolecularFeatures) -> np.ndarray:
        """
        Prepare features using explicit feature mapping (recommended method).
        
        Uses the feature mapping to place each feature at its exact expected
        position, ensuring consistency across different environments and
        versions. Handles missing features gracefully with zero-filling.
        
        Args:
            mol_features: MolecularFeatures object with named features
            
        Returns:
            np.ndarray: Feature array with exact dimension matching training
            
        Process:
            1. Create feature name->value mapping from molecular features
            2. Allocate zero-filled array of correct dimension
            3. Place each feature at its mapped index
            4. Apply preprocessing transformations
        """
        target_registry = self.feature_mapping['features']
        input_dim = self.feature_mapping['metadata']['input_dimension']
        
        features_dict: Dict[str, float] = {}
        
        # Process descriptors with names
        if mol_features.descriptors is not None and mol_features.descriptor_names:
            descriptors = ensure_numpy_array(mol_features.descriptors)
            
            for i, name in enumerate(mol_features.descriptor_names):
                if i < len(descriptors):
                    features_dict[name] = descriptors[i]
        
        # Process each fingerprint type
        fp_types = [
            ('morgan_fingerprint', 'morgan_fingerprint_'),
            ('morgan_feature_fp', 'morgan_feature_fp_'),
            ('rdkit_fingerprint', 'rdkit_fingerprint_'),
            ('avalon_fingerprint', 'avalon_fingerprint_'),
            ('pattern_fingerprint', 'pattern_fingerprint_'),
            ('layered_fingerprint', 'layered_fingerprint_')
        ]
        
        for fp_key, fp_prefix in fp_types:
            if fp_key in mol_features.fingerprints and mol_features.fingerprints[fp_key] is not None:
                fp_data = ensure_numpy_array(mol_features.fingerprints[fp_key])
                for i in range(len(fp_data)):
                    features_dict[f"{fp_prefix}{i}"] = fp_data[i]
        
        # Allocate the raw feature vector
        X_raw = np.zeros((1, input_dim), dtype=np.float32)
        
        # Place every available value at its original index
        for feat_name, val in features_dict.items():
            info = target_registry.get(feat_name)
            if info is not None:
                idx = info['original_index']
                X_raw[0, idx] = val
        
        # Run the saved pre-processor
        if self.preprocessor:
            X_proc = self.preprocessor.transform(X_raw)
        else:
            X_proc = X_raw
        
        return X_proc 