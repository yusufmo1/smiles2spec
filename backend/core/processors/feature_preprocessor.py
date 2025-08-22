"""Feature preprocessing and scaling for machine learning."""
import numpy as np
from sklearn.preprocessing import StandardScaler
import json
from joblib import load, dump
from typing import Optional

from ...config import settings
from ...utils import logger

class FeaturePreprocessor:
    """
    Feature preprocessing and scaling for machine learning models.
    
    Handles feature transformation, scaling, and cleaning to prepare
    molecular features for mass spectrum prediction. Includes automatic
    detection and handling of binary features, missing values, and
    low-variance features.
    
    Key operations:
    - Sign-preserving log transform for continuous features
    - Automatic binary feature detection
    - NaN and infinity handling with multiple strategies
    - Variance-based feature filtering
    - Separate scaling for binary and continuous features
    
    Attributes:
        config: Processing configuration dictionary
        binary_feature_mask: Boolean mask for binary features
        valid_feature_mask: Boolean mask for valid features after filtering
        scaler: StandardScaler for continuous features
        binary_scaler: Optional StandardScaler for binary features
        nan_mask: Mask for features without NaN/inf values
    """
    
    def __init__(self, config=None):
        self.config = config if config is not None else settings.feature_processing.__dict__
        self.binary_feature_mask = None
        self.valid_feature_mask = None
        self.scaler = None
        self.binary_scaler = None
        self.nan_mask = None
    
    def _sign_preserving_log_transform(self, X):
        """
        Apply sign-preserving log transform: sign(x)*log1p(|x|).
        
        This transformation helps normalize skewed features while preserving
        sign information, which is important for features like partial charges.
        
        Args:
            X: Feature array
            
        Returns:
            Transformed array with same shape and sign structure
        """
        return np.sign(X) * np.log1p(np.abs(X))
    
    def fit(self, X):
        """
        Fit the preprocessor to training data.
        
        Learns feature statistics, identifies binary features, and
        prepares scalers for transformation. Should be called once
        on training data before transform().
        
        Args:
            X: Training feature matrix (n_samples, n_features)
            
        Returns:
            self: Fitted preprocessor instance
            
        Process:
            1. Handle NaN/inf values
            2. Detect binary features
            3. Create valid feature mask
            4. Fit scalers on valid features
        """
        if X is None or len(X) == 0:
            logger.warning("Empty feature matrix, skipping preprocessing")
            return self
        X = np.array(X) if not isinstance(X, np.ndarray) else X
        
        self._handle_nans(X)
        if self.config.get('auto_detect_binary', True):
            self._detect_binary_features(X)
        self._create_valid_feature_mask(X)
        
        if self.valid_feature_mask is not None:
            X_valid = X[:, self.valid_feature_mask]
            self._scale_features(X_valid, fit=True)
        
        return self
    
    def transform(self, X):
        """
        Transform new data using fitted preprocessor.
        
        Applies learned transformations to new feature data. Must be
        called after fit() to ensure proper preprocessing.
        
        Args:
            X: Feature matrix to transform (n_samples, n_features)
            
        Returns:
            Transformed and scaled feature matrix
            
        Transformations applied:
            1. Handle NaN/inf values
            2. Log transform continuous features
            3. Select valid features
            4. Apply scaling
        """
        if X is None or len(X) == 0:
            logger.warning("Empty feature matrix, skipping preprocessing")
            return X
        
        X = np.array(X) if not isinstance(X, np.ndarray) else X
        X = self._handle_nan_values(X)
        
        # Apply the sign-preserving log transform to continuous features.
        if self.binary_feature_mask is not None:
            continuous_indices = np.where(~self.binary_feature_mask)[0]
            if len(continuous_indices) > 0:
                X[:, continuous_indices] = self._sign_preserving_log_transform(X[:, continuous_indices])
        
        if self.valid_feature_mask is None:
            logger.error("FeaturePreprocessor not fitted before transform")
            return X
        
        X_valid = X[:, self.valid_feature_mask]
        X_scaled = self._scale_features(X_valid, fit=False)
        return X_scaled
    
    def _handle_nans(self, X):
        """
        Check for and handle NaN and infinity values during fitting.
        
        Identifies features with invalid values and creates masks for
        feature selection. Logs statistics about invalid values found.
        
        Args:
            X: Feature matrix to check
            
        Returns:
            Cleaned feature matrix
            
        Side effects:
            Sets self.nan_mask for features without invalid values
        """
        # Check for NaN values
        nan_counts = np.isnan(X).sum(axis=0)
        # Check for infinity values
        inf_counts = np.isinf(X).sum(axis=0)
        
        total_invalid = nan_counts + inf_counts
        if total_invalid.sum() > 0:
            logger.info(f"Found {nan_counts.sum()} NaN values and {inf_counts.sum()} infinity values across {(total_invalid > 0).sum()} features")
        
        # Create mask for valid values (no NaN or inf)
        self.nan_mask = total_invalid == 0
        
        # Handle invalid values
        X = self._handle_nan_values(X)
        return X
    
    def _handle_nan_values(self, X):
        """
        Handle both NaN and infinity values based on configured strategy.
        
        Supports multiple strategies for dealing with invalid values:
        - 'fill_zero': Replace with zeros
        - 'fill_mean': Replace with feature means
        - 'drop_feature': Features with NaN/inf are filtered out
        
        Args:
            X: Feature matrix with potential invalid values
            
        Returns:
            Feature matrix with invalid values handled
        """
        strategy = self.config.get('handle_nan_strategy', 'drop_feature')
        X_copy = X.copy()
        if strategy == 'fill_zero':
            # Handle both NaN and infinity
            X_copy = np.nan_to_num(X_copy, nan=0.0, posinf=0.0, neginf=0.0)
        elif strategy == 'fill_mean':
            # Calculate means ignoring both NaN and infinity
            finite_mask = np.isfinite(X_copy)
            col_means = np.zeros(X_copy.shape[1])
            
            for col_idx in range(X_copy.shape[1]):
                col_data = X_copy[:, col_idx]
                finite_col_data = col_data[finite_mask[:, col_idx]]
                if len(finite_col_data) > 0:
                    col_means[col_idx] = np.mean(finite_col_data)
            
            # Replace both NaN and infinity with means
            non_finite_mask = ~np.isfinite(X_copy)
            for col_idx in range(X_copy.shape[1]):
                X_copy[non_finite_mask[:, col_idx], col_idx] = col_means[col_idx]
        
        return X_copy
    
    def _detect_binary_features(self, X):
        """
        Automatically detect binary features in the dataset.
        
        Identifies features that only contain values in {0, 1} or subsets.
        Binary features are handled differently during scaling to preserve
        their discrete nature.
        
        Args:
            X: Feature matrix to analyze
            
        Side effects:
            Sets self.binary_feature_mask
        """
        unique_vals = [np.unique(X[:, i][np.isfinite(X[:, i])]) for i in range(X.shape[1])]
        is_binary = [len(u) <= 2 and set(u).issubset({0, 1}) for u in unique_vals]
        self.binary_feature_mask = np.array(is_binary)
        binary_count = np.sum(self.binary_feature_mask)
        if binary_count > 0:
            logger.info(f"Detected {binary_count} binary features")
    
    def _create_valid_feature_mask(self, X):
        """
        Create mask for valid features based on multiple criteria.
        
        Filters out features based on:
        - NaN/inf values (if drop_feature strategy)
        - Low variance threshold
        
        Args:
            X: Feature matrix to analyze
            
        Returns:
            Boolean mask of valid features
            
        Side effects:
            Sets self.valid_feature_mask
        """
        valid_mask = np.ones(X.shape[1], dtype=bool)
        strategy = self.config.get('handle_nan_strategy', 'drop_feature')
        if strategy == 'drop_feature' and self.nan_mask is not None:
            valid_mask = valid_mask & self.nan_mask
            nan_count = np.sum(np.isnan(X).any(axis=0))
            inf_count = np.sum(np.isinf(X).any(axis=0))
            logger.info(f"Dropped {np.sum(~self.nan_mask)} features with invalid values ({nan_count} with NaN, {inf_count} with infinity)")
        
        variance_threshold = self.config.get('min_variance_threshold', 1e-8)
        if variance_threshold > 0:
            # Replace infinity with NaN so nanvar can handle it
            X_for_var = X.copy()
            X_for_var[np.isinf(X_for_var)] = np.nan
            variances = np.nanvar(X_for_var, axis=0)
            variances = np.nan_to_num(variances, nan=0.0)  # Replace any NaN variance with 0
            
            high_variance_mask = variances > variance_threshold
            valid_mask = valid_mask & high_variance_mask
            logger.info(f"Dropped {np.sum(~high_variance_mask)} features with variance below {variance_threshold}")
            
            # Free memory
            del X_for_var
            del variances
        
        self.valid_feature_mask = valid_mask
        logger.info(f"Final feature count: {np.sum(valid_mask)} of {X.shape[1]} original features")
        return valid_mask
    
    def _scale_features(self, X, fit=True):
        """
        Scale features using appropriate methods for each feature type.
        
        Applies StandardScaler to continuous features and optionally to
        binary features based on configuration. Handles edge cases and
        ensures numerical stability.
        
        Args:
            X: Valid features to scale
            fit: Whether to fit scalers (True) or just transform (False)
            
        Returns:
            Scaled feature matrix as float32 for memory efficiency
            
        Note:
            Binary and continuous features can be scaled separately
        """
        if X is None or X.shape[0] == 0 or X.shape[1] == 0:
            return X
        scale_continuous = self.config.get('scaling', {}).get('scale_continuous', True)
        scale_binary = self.config.get('scaling', {}).get('scale_binary', False)
        
        # If no scaling is required, clip and return.
        if not scale_continuous and not scale_binary:
            X_clipped = np.clip(X, np.finfo(np.float32).min, np.finfo(np.float32).max)
            return X_clipped.astype(np.float32)
        
        # If binary feature detection is off or not available.
        if self.binary_feature_mask is None or not self.config.get('auto_detect_binary', True):
            if scale_continuous:
                if fit:
                    self.scaler = StandardScaler()
                    X_scaled = self.scaler.fit_transform(X)
                else:
                    if self.scaler is None:
                        logger.warning("Scaler not fitted before transform. Returning unscaled features.")
                        X_scaled = X
                    else:
                        X_scaled = self.scaler.transform(X)
            else:
                X_scaled = X
            X_clipped = np.clip(X_scaled, np.finfo(np.float32).min, np.finfo(np.float32).max)
            return X_clipped.astype(np.float32)
        
        # When binary feature mask is available.
        binary_valid = self.binary_feature_mask[self.valid_feature_mask]
        X_processed = X.copy()
        if scale_continuous:
            continuous_mask = ~binary_valid
            if np.any(continuous_mask):
                continuous_indices = np.where(continuous_mask)[0]
                continuous_data = X[:, continuous_indices]
                if fit:
                    self.scaler = StandardScaler()
                    X_processed[:, continuous_indices] = self.scaler.fit_transform(continuous_data)
                else:
                    if self.scaler is None:
                        logger.warning("Scaler not fitted before transform. Returning unscaled continuous features.")
                    else:
                        X_processed[:, continuous_indices] = self.scaler.transform(continuous_data)
                del continuous_data
        if scale_binary and np.any(binary_valid):
            binary_indices = np.where(binary_valid)[0]
            if fit:
                self.binary_scaler = StandardScaler()
                X_processed[:, binary_indices] = self.binary_scaler.fit_transform(X[:, binary_indices])
            else:
                if hasattr(self, 'binary_scaler') and self.binary_scaler is not None:
                    X_processed[:, binary_indices] = self.binary_scaler.transform(X[:, binary_indices])
        X_clipped = np.clip(X_processed, np.finfo(np.float32).min, np.finfo(np.float32).max)
        return X_clipped.astype(np.float32)
    
    @classmethod
    def load(cls, file_path):
        """
        Load a saved preprocessor from file.
        
        Restores all preprocessor state including fitted scalers,
        feature masks, and configuration.
        
        Args:
            file_path: Path to saved preprocessor (.joblib file)
            
        Returns:
            Loaded FeaturePreprocessor instance or None on error
        """
        try:
            loaded_state = load(file_path)
            preprocessor = cls(loaded_state.get('config', {}))
            preprocessor.binary_feature_mask = loaded_state.get('binary_feature_mask')
            preprocessor.valid_feature_mask = loaded_state.get('valid_feature_mask')
            preprocessor.scaler = loaded_state.get('scaler')
            preprocessor.binary_scaler = loaded_state.get('binary_scaler')
            preprocessor.nan_mask = loaded_state.get('nan_mask')
            preprocessor.nan_feature_indices = loaded_state.get('nan_feature_indices', [])
            preprocessor.low_variance_feature_indices = loaded_state.get('low_variance_feature_indices', [])
            logger.info(f"Loaded feature preprocessor from {file_path}")
            return preprocessor
        except Exception as e:
            logger.error(f"Error loading preprocessor from {file_path}: {str(e)}")
            return None
        
    def save(self, file_path):
        """
        Save the preprocessor state to file.
        
        Serializes all preprocessor components including scalers,
        masks, and configuration for later reuse.
        
        Args:
            file_path: Destination path for saved preprocessor
            
        Returns:
            bool: True if save successful, False otherwise
        """
        try:
            state = {
                'binary_feature_mask': self.binary_feature_mask,
                'valid_feature_mask': self.valid_feature_mask,
                'scaler': self.scaler,
                'binary_scaler': getattr(self, 'binary_scaler', None),
                'nan_mask': getattr(self, 'nan_mask', None),
                'nan_feature_indices': getattr(self, 'nan_feature_indices', []),
                'low_variance_feature_indices': getattr(self, 'low_variance_feature_indices', []),
                'config': self.config
            }
            dump(state, file_path)
            logger.info(f"Saved feature preprocessor to {file_path}")
            return True
        except Exception as e:
            logger.error(f"Error saving preprocessor to {file_path}: {str(e)}")
            return False 