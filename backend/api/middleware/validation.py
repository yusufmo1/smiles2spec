"""Validation middleware."""
from functools import wraps
from flask import request, jsonify
from pydantic import BaseModel, ValidationError
from ...utils import logger

def validate_json(schema_class: BaseModel):
    """Decorator for JSON validation using Pydantic schemas."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            try:
                data = request.get_json()
                if data is None:
                    return jsonify({"error": "No JSON data provided"}), 400
                
                # Debug logging
                logger.debug(f"Validating data with {schema_class.__name__}: {data}")
                
                validated_data = schema_class(**data)
                return f(validated_data, *args, **kwargs)
            except ValidationError as e:
                logger.warning(f"Validation error for {schema_class.__name__}: {e}")
                logger.debug(f"Raw data that caused validation error: {data}")
                return jsonify({
                    "error": "Validation failed", 
                    "details": e.errors(),
                    "schema": schema_class.__name__
                }), 400
            except Exception as e:
                logger.error(f"Validation middleware error: {str(e)}")
                return jsonify({"error": "Internal server error"}), 500
        return decorated_function
    return decorator

def validate_query_params(schema_class: BaseModel):
    """Decorator for query parameter validation."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            try:
                # Convert query args to dict
                query_data = request.args.to_dict()
                validated_data = schema_class(**query_data)
                return f(validated_data, *args, **kwargs)
            except ValidationError as e:
                logger.warning(f"Query validation error: {e}")
                return jsonify({
                    "error": "Invalid query parameters", 
                    "details": e.errors()
                }), 400
            except Exception as e:
                logger.error(f"Query validation middleware error: {str(e)}")
                return jsonify({"error": "Internal server error"}), 500
        return decorated_function
    return decorator 