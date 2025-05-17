"""
Model loading and prediction functionality.
"""

import os
import numpy as np
import json
from typing import Dict
from joblib import load
from utils import logger, handle_error, ensure_numpy_array, convert_np_to_list

class ModelHandler:
    def __init__(self):
        self.model = None
        self.preprocessor = None
        self.feature_mapping = None
    
    def load_model(self, model_path):
        """Load a trained model from file."""
        try:
            if os.path.exists(model_path):
                self.model = load(model_path)
                logger.info(f"Loaded model from {model_path}")
                return True
            else:
                logger.error(f"Model file not found: {model_path}")
                return False
        except Exception as e:
            error_result = handle_error(e, f"loading model from {model_path}")
            logger.error(f"Error loading model: {error_result['error']}")
            return False
    
    def load_preprocessor(self, preprocessor_path):
        """Load a feature preprocessor from file."""
        from feature_preprocessor import FeaturePreprocessor
        
        try:
            if os.path.exists(preprocessor_path):
                self.preprocessor = FeaturePreprocessor.load(preprocessor_path)
                logger.info(f"Loaded preprocessor from {preprocessor_path}")
                return True
            else:
                logger.error(f"Preprocessor file not found: {preprocessor_path}")
                return False
        except Exception as e:
            error_result = handle_error(e, f"loading preprocessor from {preprocessor_path}")
            logger.error(f"Error loading preprocessor: {error_result['error']}")
            return False
    
    def load_feature_mapping(self, mapping_path):
        """Load feature mapping from file."""
        try:
            if os.path.exists(mapping_path):
                with open(mapping_path, 'r') as f:
                    self.feature_mapping = json.load(f)
                logger.info(f"Loaded feature mapping from {mapping_path}")
                return True
            else:
                logger.error(f"Feature mapping file not found: {mapping_path}")
                return False
        except Exception as e:
            error_result = handle_error(e, f"loading feature mapping from {mapping_path}")
            logger.error(f"Error loading feature mapping: {error_result['error']}")
            return False
    
    def initialize(self, model_path, preprocessor_path, mapping_path=None):
        """Initialize the model handler with all required components."""
        model_loaded = self.load_model(model_path)
        preprocessor_loaded = self.load_preprocessor(preprocessor_path)
        
        if mapping_path:
            mapping_loaded = self.load_feature_mapping(mapping_path)
        else:
            mapping_loaded = True  # Skip if not provided
        
        return model_loaded and preprocessor_loaded and mapping_loaded
    
    def prepare_features(self, mol_features):
        """Extract and preprocess features from molecular features dictionary."""
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
            error_result = handle_error(e, "preparing features for prediction")
            raise ValueError(f"Error preparing features: {error_result['error']}")
    
    def _prepare_features_without_mapping(self, mol_features):
        """Original feature preparation method without mapping."""
        # Extract features in the right order
        features = []
        if 'descriptors' in mol_features and mol_features['descriptors'] is not None:
            features.append(ensure_numpy_array(mol_features['descriptors']))
        if 'morgan_fingerprint' in mol_features and mol_features['morgan_fingerprint'] is not None:
            features.append(ensure_numpy_array(mol_features['morgan_fingerprint']))
        
        additional_keys = ['morgan_feature_fp', 'atom_pairs', 'electronic_features', 'substructure_features',
                        'rdkit_fingerprint', 'avalon_fingerprint', 'pattern_fingerprint', 'layered_fingerprint']
        for key in additional_keys:
            if key in mol_features and mol_features[key] is not None:
                features.append(ensure_numpy_array(mol_features[key]))
        
        if not features:
            raise ValueError("No valid features found in molecular data")
        
        # Concatenate features
        X = np.concatenate(features).reshape(1, -1)
        
        # Preprocess if preprocessor is available
        if self.preprocessor:
            X = self.preprocessor.transform(X)
        
        return X
    
    def _prepare_features_with_mapping(self, mol_features):
        """Prepare features using the feature mapping for consistent ordering."""
        # ------------------------------------------------------------------
        # 1) Build a dict name -> value for *all* possible raw features
        # ------------------------------------------------------------------
        target_registry = self.feature_mapping['features']   # full registry
        input_dim = self.feature_mapping['metadata']['input_dimension']  # 9 416

        features_dict: Dict[str, float] = {}
        
        # Process descriptors with names
        if 'descriptors' in mol_features and 'descriptor_names' in mol_features:
            descriptors = ensure_numpy_array(mol_features['descriptors'])
            descriptor_names = mol_features['descriptor_names']
            
            for i, name in enumerate(descriptor_names):
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
            if fp_key in mol_features and mol_features[fp_key] is not None:
                fp_data = ensure_numpy_array(mol_features[fp_key])
                for i in range(len(fp_data)):
                    features_dict[f"{fp_prefix}{i}"] = fp_data[i]
        
        # ------------------------------------------------------------------
        # 2) Allocate the *raw* feature vector (9 416 columns)
        # ------------------------------------------------------------------
        X_raw = np.zeros((1, input_dim), dtype=np.float32)

        # 3) Place every available value at its ORIGINAL index
        for feat_name, val in features_dict.items():
            info = target_registry.get(feat_name)
            if info is not None:
                idx = info['original_index']
                X_raw[0, idx] = val

        # 4) Run the saved pre-processor (drops → log → scaling)
        if self.preprocessor:
            X_proc = self.preprocessor.transform(X_raw)
        else:
            X_proc = X_raw

        return X_proc
    
    def predict(self, mol_features):
        """Make a prediction using the loaded model."""
        if self.model is None:
            raise ValueError("Model not loaded. Call initialize() first.")
        
        try:
            # Prepare features
            X = self.prepare_features(mol_features)
            
            # Make prediction
            prediction = self.model.predict(X)
            
            # Return as 1D array
            return prediction[0]
        except Exception as e:
            error_result = handle_error(e, "making prediction")
            raise ValueError(f"Error making prediction: {error_result['error']}")