"""
Centralized Application Settings

Configuration management for the SMILES2SPEC backend application.
Uses dataclasses for type-safe configuration with defaults and
environment variable overrides.

Environment Variables:
- API_HOST: Server host (default: 0.0.0.0)
- API_PORT: Server port (default: 5050)
- API_DEBUG: Debug mode (default: True)
- OPENROUTER_API_KEY: API key for LLM services

Configuration is loaded from:
1. Default values in dataclasses
2. .env file in backend directory
3. Environment variable overrides
"""
import os
from typing import Dict, Any
from dataclasses import dataclass, field
from dotenv import load_dotenv

# Get the backend directory
BACKEND_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Load environment variables from .env file
dotenv_path = os.path.join(BACKEND_DIR, '.env')
if os.path.exists(dotenv_path):
    load_dotenv(dotenv_path)
    print(f"Loaded environment variables from {dotenv_path}")
else:
    print(f"No .env file found at {dotenv_path}")

@dataclass
class DatabaseConfig:
    """
    Database configuration (placeholder for future use).
    Currently unused as the application doesn't require a database.
    """
    pass

@dataclass
class APIConfig:
    """
    API server configuration.
    
    Attributes:
        host: Server bind address
        port: Server port number
        debug: Enable debug mode with hot reload
        model_path: Path to trained ML model file
        preprocessor_path: Path to feature preprocessor
        feature_mapping_path: Path to feature index mapping
    """
    host: str = '0.0.0.0'
    port: int = 5050
    debug: bool = True
    model_path: str = field(default_factory=lambda: os.path.join(BACKEND_DIR, 'models/spectrum_predictor.pkl'))
    preprocessor_path: str = field(default_factory=lambda: os.path.join(BACKEND_DIR, 'models/feature_preprocessor.pkl'))
    feature_mapping_path: str = field(default_factory=lambda: os.path.join(BACKEND_DIR, 'models/feature_mapping.json'))

@dataclass
class FeatureConfig:
    """
    Molecular feature extraction configuration.
    
    Controls which molecular features are extracted from SMILES:
    - RDKit descriptors (200+ properties)
    - Multiple fingerprint types with configurable sizes
    - Atom/bond counting statistics
    - Electronic features (partial charges, etc.)
    """
    extract_descriptors: bool = True
    fingerprints: Dict[str, Any] = field(default_factory=lambda: {
        'morgan': {'enabled': True, 'radii': [1, 2, 3], 'size': 1024},
        'morgan_feature': {'enabled': True, 'radius': 2, 'size': 1024},
        'maccs': {'enabled': True},
        'topological': {'enabled': True, 'size': 1024},
        'rdkit': {'enabled': True, 'size': 2048},
        'avalon': {'enabled': True, 'size': 1024},
        'pattern': {'enabled': True, 'size': 1024},
        'layered': {'enabled': True, 'size': 2048}
    })
    count_features: Dict[str, bool] = field(default_factory=lambda: {
        'bond_counts': True,
        'atom_counts': True
    })
    extract_electronic: bool = True
    extract_substructures: bool = False

@dataclass
class FeatureProcessingConfig:
    """
    Feature preprocessing and scaling configuration.
    
    Attributes:
        handle_nan_strategy: How to handle missing values
            - 'drop_feature': Remove features with NaN
            - 'fill_zero': Replace with zeros
            - 'fill_mean': Replace with feature mean
        min_variance_threshold: Remove features below this variance
        auto_detect_binary: Automatically detect binary features
        scaling: Separate scaling options for binary/continuous
    """
    handle_nan_strategy: str = 'drop_feature'
    min_variance_threshold: float = 1e-8
    auto_detect_binary: bool = True
    scaling: Dict[str, bool] = field(default_factory=lambda: {
        'scale_continuous': True,
        'scale_binary': False
    })

@dataclass
class SpectralConfig:
    """
    Mass spectrum processing configuration.
    
    Attributes:
        max_peaks: Maximum number of peaks to retain
        bin_size: m/z binning resolution (Da)
        max_mz: Maximum m/z value to consider
        sort_peaks_by: Peak sorting criterion ('intensity' or 'mz')
        intensity_distribution_bins: Bins for intensity histogram
        peak_distribution_bins: Bins for m/z distribution
    """
    max_peaks: int = 499
    bin_size: float = 1.0
    max_mz: int = 499
    sort_peaks_by: str = 'intensity'
    intensity_distribution_bins: int = 20
    peak_distribution_bins: int = 50

@dataclass
class LLMConfig:
    """
    Large Language Model integration configuration.
    
    Configuration for OpenRouter API which provides access to
    various LLMs for chat and SMILES generation features.
    
    Attributes:
        openrouter_api_key: API key (set via environment)
        openrouter_base_url: API endpoint
        site_url: Your application URL for tracking
        site_name: Application name for OpenRouter dashboard
        models: Available models for different use cases
    """
    openrouter_api_key: str = field(default_factory=lambda: os.getenv('OPENROUTER_API_KEY', ''))
    openrouter_base_url: str = "https://openrouter.ai/api/v1"
    site_url: str = "https://smiles2spec.app"
    site_name: str = "SMILES2Spec App"
    models: Dict[str, str] = field(default_factory=lambda: {
        'default': 'google/gemma-3-27b-it:free',
        'fast': 'google/gemma-3-27b-it:free',
        'smart': 'google/gemma-3-27b-it:freet'
    })

@dataclass
class Settings:
    """
    Main application settings container.
    
    Aggregates all configuration sections and provides
    a unified interface for accessing settings throughout
    the application.
    """
    api: APIConfig = field(default_factory=APIConfig)
    features: FeatureConfig = field(default_factory=FeatureConfig)
    feature_processing: FeatureProcessingConfig = field(default_factory=FeatureProcessingConfig)
    spectral: SpectralConfig = field(default_factory=SpectralConfig)
    llm: LLMConfig = field(default_factory=LLMConfig)
    
    @classmethod
    def load(cls) -> 'Settings':
        """
        Load settings from environment and config files.
        
        Loading order (later overrides earlier):
        1. Default values from dataclass definitions
        2. Values from .env file if present
        3. Environment variable overrides
        
        Returns:
            Settings instance with loaded configuration
        """
        settings = cls()
        
        # Override with environment variables if present
        if os.getenv('API_HOST'):
            settings.api.host = os.getenv('API_HOST')
        if os.getenv('API_PORT'):
            settings.api.port = int(os.getenv('API_PORT'))
        if os.getenv('API_DEBUG'):
            settings.api.debug = os.getenv('API_DEBUG').lower() == 'true'
        
        # Load LLM configuration from environment
        if os.getenv('OPENROUTER_API_KEY'):
            settings.llm.openrouter_api_key = os.getenv('OPENROUTER_API_KEY')
            
        return settings

# Global settings instance - import this in other modules
# Example: from config.settings import settings
settings = Settings.load() 