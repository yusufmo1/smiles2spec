"""API routes package."""
from .prediction import prediction_bp
from .export import export_bp
from .chat import chat_bp
from .health import health_bp
from .upload import upload_bp

__all__ = [
    'prediction_bp',
    'export_bp', 
    'chat_bp',
    'health_bp',
    'upload_bp'
] 