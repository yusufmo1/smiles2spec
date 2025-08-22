"""Export routes."""
import time
from flask import Blueprint, request, make_response, jsonify
from ...services.prediction_service import PredictionService
from ...utils.formats import peaks_to_msp
from ...utils import logger
from ..schemas.prediction import ExportMSPRequest, ExportMSPBatchRequest
from ..middleware.validation import validate_json

export_bp = Blueprint('export', __name__, url_prefix='')
prediction_service = PredictionService()

@export_bp.route('/export_msp', methods=['POST'])
@validate_json(ExportMSPRequest)
def export_msp(validated_data: ExportMSPRequest):
    """
    Export a single predicted mass spectrum as MSP file.
    
    The MSP (Mass Spectral Peak) format is a standard text format for mass spectrometry
    data exchange, commonly used in spectral libraries and databases.
    
    Args:
        validated_data: ExportMSPRequest containing:
            - smiles: SMILES string of the molecule to predict and export
            
    Returns:
        Text file response with:
            - Content-Type: text/plain
            - MSP formatted spectral data
            - Automatic download with descriptive filename
            
    MSP Format includes:
        - NAME: Molecule identifier
        - MW: Molecular weight
        - FORMULA: Molecular formula
        - SMILES: Input SMILES string
        - Num Peaks: Number of spectral peaks
        - Peak list: m/z intensity pairs
        
    Raises:
        500: Error during prediction or MSP generation
    """
    try:
        result = prediction_service.predict_spectrum(validated_data.smiles)
        msp_txt, fname = peaks_to_msp(result.to_dict())
        
        resp = make_response(msp_txt)
        resp.headers["Content-Type"] = "text/plain; charset=utf-8"
        resp.headers["Content-Disposition"] = f'attachment; filename="{fname}.msp"'
        return resp
    except Exception as e:
        logger.error(f"MSP export error: {str(e)}")
        return jsonify({"error": str(e)}), 500

@export_bp.route('/export_msp_batch', methods=['POST'])
@validate_json(ExportMSPBatchRequest)
def export_msp_batch(validated_data: ExportMSPBatchRequest):
    """
    Export multiple predicted mass spectra as a single MSP file.
    
    Processes a batch of SMILES strings and combines all successfully predicted
    spectra into a single MSP file for bulk analysis or library creation.
    
    Args:
        validated_data: ExportMSPBatchRequest containing:
            - smiles_list: List of SMILES strings to process
            
    Returns:
        Text file response with:
            - Combined MSP data for all valid molecules
            - Timestamped filename for unique identification
            - Automatic download trigger
            
    Notes:
        - Invalid SMILES are skipped with error logging
        - At least one valid spectrum required for successful export
        - Each spectrum is separated by double newlines in the output
        
    Raises:
        400: No valid spectra could be generated
        500: Error during batch processing
    """
    try:
        # Generate MSP for each SMILES
        msp_entries = []
        timestamp = int(time.time())
        
        for smiles in validated_data.smiles_list:
            try:
                result = prediction_service.predict_spectrum(smiles)
                msp_txt, _ = peaks_to_msp(result.to_dict())
                msp_entries.append(msp_txt)
            except Exception as e:
                logger.error(f"Error processing SMILES {smiles}: {str(e)}")
                continue
        
        if not msp_entries:
            return jsonify({"error": "No valid spectra generated"}), 400
        
        # Combine all entries
        full_msp = "\n\n".join(msp_entries)
        
        resp = make_response(full_msp)
        resp.headers["Content-Type"] = "text/plain; charset=utf-8"
        resp.headers["Content-Disposition"] = f'attachment; filename="batch_export_{timestamp}.msp"'
        return resp
    except Exception as e:
        logger.error(f"Batch MSP export error: {str(e)}")
        return jsonify({"error": str(e)}), 500 