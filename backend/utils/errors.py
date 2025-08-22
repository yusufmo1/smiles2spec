"""Centralized error handling."""
from typing import Dict, Any
import logging

logger = logging.getLogger(__name__)

class APIException(Exception):
    """Base API exception."""
    def __init__(self, message: str, status_code: int = 400, details: Dict[str, Any] = None):
        self.message = message
        self.status_code = status_code
        self.details = details or {}
        super().__init__(self.message)

class ValidationError(APIException):
    """Input validation error."""
    def __init__(self, message: str, details: Dict[str, Any] = None):
        super().__init__(message, 400, details)

class SMILESError(APIException):
    """SMILES processing error."""
    def __init__(self, message: str, smiles: str = None):
        details = {"smiles": smiles} if smiles else {}
        super().__init__(message, 400, details)

class ModelError(APIException):
    """ML model error."""
    def __init__(self, message: str, details: Dict[str, Any] = None):
        super().__init__(message, 500, details)

class ExternalServiceError(APIException):
    """External service error."""
    def __init__(self, message: str, service: str = None):
        details = {"service": service} if service else {}
        super().__init__(message, 503, details)

class FileProcessingError(APIException):
    """File processing error."""
    def __init__(self, message: str, filename: str = None):
        details = {"filename": filename} if filename else {}
        super().__init__(message, 400, details) 