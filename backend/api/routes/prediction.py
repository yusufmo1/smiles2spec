"""Prediction routes."""
from flask import Blueprint, jsonify
from ...services.prediction_service import PredictionService
from ...utils import logger
from ..schemas.prediction import PredictRequest
from ..middleware.validation import validate_json

prediction_bp = Blueprint('prediction', __name__, url_prefix='')
prediction_service = PredictionService()

@prediction_bp.route('/predict', methods=['POST'])
@validate_json(PredictRequest)
def predict(validated_data: PredictRequest):
    """
    Predict mass spectrum from SMILES molecular structure.
    
    This endpoint accepts a SMILES string and returns predicted mass spectrum data
    including m/z values, intensities, and molecular properties.
    
    Args:
        validated_data: PredictRequest object containing the SMILES string
        
    Returns:
        JSON response containing:
        - mz: Array of mass-to-charge ratios
        - intensity: Array of intensity values (0-100)
        - smiles: Input SMILES string
        - molecular_weight: Calculated molecular weight
        - molecular_formula: Chemical formula
        - Additional molecular properties
        
    Raises:
        400: Invalid SMILES string
        500: Internal server error during prediction
    """
    try:
        result = prediction_service.predict_spectrum(validated_data.smiles)
        return jsonify(result.to_dict())
    except Exception as e:
        logger.error(f"Prediction error: {str(e)}")
        return jsonify({"error": str(e)}), 500

 