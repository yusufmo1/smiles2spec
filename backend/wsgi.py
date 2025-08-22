"""
WSGI Entry Point for SMILES2SPEC Backend

This module provides the WSGI application object for production deployment
with Gunicorn. It imports the Flask application from the main module and
exposes it for WSGI servers.

Usage:
    gunicorn --config backend/gunicorn.conf.py backend.wsgi:app

The application includes:
- Flask API with all routes
- CORS configuration for production
- ML model initialization
- Error handling middleware
"""

# Import the Flask application from the backend module
from backend.main import create_app
from backend.utils import logger

# Create the WSGI application
app = create_app()

# Log successful initialization
logger.info("WSGI application initialized successfully")

if __name__ == "__main__":
    # This won't be called in production, but useful for testing
    logger.warning("Running WSGI app directly - use Gunicorn for production")
    app.run(host="0.0.0.0", port=5050, debug=False)