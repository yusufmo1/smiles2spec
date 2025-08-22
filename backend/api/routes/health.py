"""Health check routes."""
from flask import Blueprint, jsonify
from ...services.prediction_service import PredictionService
from ...utils import logger

health_bp = Blueprint('health', __name__, url_prefix='')

@health_bp.route('/health', methods=['GET'])
def health_check():
    """
    Service health check endpoint.
    
    Provides comprehensive health status information about the API service,
    including model readiness and service availability. Used for monitoring
    and load balancer health checks.
    
    Returns:
        JSON response containing:
            - status: Overall health status ('healthy', 'degraded', 'unhealthy')
            - service_ready: Whether the prediction service is operational
            - version: API version number
            - architecture: Service architecture type
            - model_loaded: Whether ML models are loaded and ready
            
    Status Codes:
        200: Service is healthy or degraded but functional
        500: Service is unhealthy and cannot process requests
        
    Notes:
        - 'degraded' status indicates service is running but models aren't loaded
        - Used by Docker health checks and monitoring systems
    """
    try:
        prediction_service = PredictionService()
        
        status = {
            "status": "healthy",
            "service_ready": prediction_service.is_ready(),
            "version": "2.0.0",
            "architecture": "service-oriented"
        }
        
        # Add more detailed health information
        if prediction_service.is_ready():
            status["model_loaded"] = True
        else:
            status["model_loaded"] = False
            status["status"] = "degraded"
        
        return jsonify(status)
    except Exception as e:
        logger.error(f"Health check error: {str(e)}")
        return jsonify({
            "status": "unhealthy",
            "error": str(e)
        }), 500 