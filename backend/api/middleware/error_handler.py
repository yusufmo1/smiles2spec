"""Global error handling middleware."""
from flask import jsonify, current_app
from ...utils.errors import APIException, ValidationError, SMILESError, ModelError, ExternalServiceError

def register_error_handlers(app):
    """Register global error handlers with the Flask app."""
    
    @app.errorhandler(APIException)
    def handle_api_exception(error):
        """Handle custom API exceptions."""
        response = {
            "error": error.message,
            "error_type": error.__class__.__name__,
            "status_code": error.status_code
        }
        if error.details:
            response["details"] = error.details
        
        current_app.logger.warning(f"API Exception: {error.message}")
        return jsonify(response), error.status_code
    
    @app.errorhandler(ValidationError)
    def handle_validation_error(error):
        """Handle validation errors."""
        response = {
            "error": error.message,
            "error_type": "ValidationError",
            "status_code": 400
        }
        if error.details:
            response["details"] = error.details
        
        current_app.logger.warning(f"Validation Error: {error.message}")
        return jsonify(response), 400
    
    @app.errorhandler(SMILESError)
    def handle_smiles_error(error):
        """Handle SMILES-related errors."""
        response = {
            "error": error.message,
            "error_type": "SMILESError",
            "status_code": 400
        }
        if error.details:
            response["details"] = error.details
        
        current_app.logger.warning(f"SMILES Error: {error.message}")
        return jsonify(response), 400
    
    @app.errorhandler(ModelError)
    def handle_model_error(error):
        """Handle model-related errors."""
        response = {
            "error": error.message,
            "error_type": "ModelError",
            "status_code": 500
        }
        if error.details:
            response["details"] = error.details
        
        current_app.logger.error(f"Model Error: {error.message}")
        return jsonify(response), 500
    
    @app.errorhandler(ExternalServiceError)
    def handle_external_service_error(error):
        """Handle external service errors."""
        response = {
            "error": error.message,
            "error_type": "ExternalServiceError", 
            "status_code": 503
        }
        if error.details:
            response["details"] = error.details
        
        current_app.logger.error(f"External Service Error: {error.message}")
        return jsonify(response), 503
    
    @app.errorhandler(404)
    def handle_not_found(error):
        """Handle 404 errors."""
        response = {
            "error": "Endpoint not found",
            "error_type": "NotFound",
            "status_code": 404
        }
        current_app.logger.warning(f"404 Error: {error}")
        return jsonify(response), 404
    
    @app.errorhandler(405)
    def handle_method_not_allowed(error):
        """Handle 405 errors."""
        response = {
            "error": "Method not allowed",
            "error_type": "MethodNotAllowed",
            "status_code": 405
        }
        current_app.logger.warning(f"405 Error: {error}")
        return jsonify(response), 405
    
    @app.errorhandler(Exception)
    def handle_generic_exception(error):
        """Handle unexpected exceptions."""
        current_app.logger.error(f"Unhandled exception: {str(error)}", exc_info=True)
        
        # Don't leak internal details in production
        if current_app.debug:
            response = {
                "error": str(error),
                "error_type": "InternalServerError",
                "status_code": 500
            }
        else:
            response = {
                "error": "Internal server error",
                "error_type": "InternalServerError",
                "status_code": 500
            }
        
        return jsonify(response), 500 