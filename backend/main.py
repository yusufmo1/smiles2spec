"""
Main Application Entry Point

Bootstraps the Flask application with all necessary configurations,
middleware, and route handlers. Uses the application factory pattern
for better testability and configuration management.

Key Features:
- CORS configuration for production and development
- Centralized error handling
- Blueprint-based route organization
- Automatic model initialization with graceful fallback

Usage:
    Development: python main.py
    Production: gunicorn main:app
"""
from flask import Flask
from flask_cors import CORS

from .config.settings import settings
from .services.prediction_service import PredictionService
from .utils import logger
from .api.middleware.error_handler import register_error_handlers

# Import all route blueprints
from .api.routes import (
    prediction_bp,
    export_bp,
    chat_bp,
    health_bp,
    upload_bp
)

def create_app() -> Flask:
    """
    Application factory function.
    
    Creates and configures the Flask application with all necessary
    components. This pattern allows for easy testing and multiple
    app instances with different configurations.
    
    Returns:
        Configured Flask application instance
        
    Configuration includes:
        - CORS for cross-origin requests
        - Error handlers for consistent error responses
        - API route blueprints
        - ML model initialization
    """
    app = Flask(__name__)
    
    # Configure app
    app.config['JSON_AS_ASCII'] = False
    
    # Configure CORS for production and development
    # Allows frontend to communicate with API across domains
    CORS(app, 
         origins=[
             "https://spectralsimulation.com",
             "https://www.spectralsimulation.com",
             "http://localhost:5173",  # Development frontend
             "http://localhost:3000",   # Alternative development port
             "http://localhost:3001"   # Docker frontend port
         ],
         methods=['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS'],
         allow_headers=['Content-Type', 'Authorization', 'Accept'],
         supports_credentials=True)
    
    # Register error handlers
    register_error_handlers(app)
    
    # Register blueprints
    app.register_blueprint(prediction_bp)
    app.register_blueprint(export_bp)
    app.register_blueprint(chat_bp)
    app.register_blueprint(health_bp)
    app.register_blueprint(upload_bp)
    
    # Initialize prediction service
    # Attempts to load ML models; continues without them if unavailable
    prediction_service = PredictionService()
    if not prediction_service.initialize():
        logger.warning("Prediction service failed to initialize - running without model")
    
    logger.info("Application initialized successfully")
    return app

# Create application instance
# This is what WSGI servers (like gunicorn) will import
app = create_app()

if __name__ == '__main__':
    # Development server - runs when script is executed directly
    # In production, use Gunicorn via: gunicorn --config gunicorn.conf.py backend.wsgi:app
    
    import os
    
    # Check if we're in production environment
    flask_env = os.getenv('FLASK_ENV', 'development').lower()
    
    if flask_env == 'production':
        logger.warning(
            "Production environment detected. "
            "Consider using Gunicorn: 'gunicorn --config gunicorn.conf.py backend.wsgi:app'"
        )
    
    logger.info(f"Starting development server in {flask_env} mode")
    app.run(
        host=settings.api.host,
        port=settings.api.port,
        debug=settings.api.debug and flask_env != 'production'
    ) 