"""
Configuration settings for the mass spectrometry prediction API.
"""

# Feature extraction configuration
FEATURE_CONFIG = {
    # Feature extraction settings
    'extract_descriptors': True,
    
    # Fingerprint configurations
    'fingerprints': {
        'morgan': {'enabled': True, 'radii': [1, 2, 3], 'size': 1024},
        'morgan_feature': {'enabled': True, 'radius': 2, 'size': 1024},
        'maccs': {'enabled': True},
        'topological': {'enabled': True, 'size': 1024},
        'rdkit': {'enabled': True, 'size': 2048},
        'avalon': {'enabled': True, 'size': 1024},
        'pattern': {'enabled': True, 'size': 1024},
        'layered': {'enabled': True, 'size': 2048}
    },
    
    # Other features
    'count_features': {
        'bond_counts': True,
        'atom_counts': True
    },
    'extract_electronic': True,
    'extract_substructures': False
}

# Spectral processing configuration
SPECTRAL_CONFIG = {
    'max_peaks': 499,
    'bin_size': 1.0,
    'max_mz': 499,
    'sort_peaks_by': 'intensity',
    'intensity_distribution_bins': 20,
    'peak_distribution_bins': 50,
}

# Feature processing configuration
FEATURE_PROCESSING_CONFIG = {
    'handle_nan_strategy': 'drop_feature',
    'min_variance_threshold': 1e-8,
    'auto_detect_binary': True,
    'scaling': {
        'scale_continuous': True,
        'scale_binary': False
    }
}

# API configuration
API_CONFIG = {
    'model_path': 'models/spectrum_predictor.pkl',
    'preprocessor_path': 'models/feature_preprocessor.pkl',
    'feature_mapping_path': 'models/feature_mapping.json',
    'debug': True,
    'host': '0.0.0.0',
    'port': 5050
}