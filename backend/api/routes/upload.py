"""File upload routes."""
import csv
from io import StringIO
from flask import Blueprint, request, jsonify
from ...utils import logger

upload_bp = Blueprint('upload', __name__, url_prefix='')

@upload_bp.route('/smiles_bulk', methods=['POST'])
def smiles_bulk():
    """
    Process bulk SMILES upload from CSV or text files.
    
    Accepts file uploads containing multiple SMILES strings and extracts them
    for batch processing. Supports both CSV and plain text formats.
    
    File Formats:
        - CSV: Uses first column for SMILES, auto-detects headers
        - TXT: One SMILES per line, handles space-separated data
        
    Request:
        - Multipart form data with 'file' field
        - Supported extensions: .csv, .txt (or no extension)
        
    Returns:
        JSON response containing:
            - smiles: List of extracted SMILES strings
            
    Notes:
        - Automatically skips header rows that don't look like SMILES
        - For space-separated files, only first token per line is used
        - Empty lines are ignored
        - UTF-8 encoding assumed (errors ignored)
        
    Raises:
        400: No file provided
        500: Error processing file
    """
    try:
        f = request.files.get("file")
        if f is None:
            return jsonify({"error": "No file provided"}), 400
        
        content = f.stream.read().decode("utf-8", errors="ignore")
        filename = f.filename.lower() if f.filename else ""
        
        if filename.endswith('.csv'):
            # Parse CSV content
            smiles_list = []
            reader = csv.reader(StringIO(content))
            for row in reader:
                if row and row[0].strip():  # Use first column for SMILES
                    smiles_list.append(row[0].strip())
            
            # Skip header row if it doesn't look like SMILES
            if smiles_list and not _looks_like_smiles(smiles_list[0]):
                smiles_list = smiles_list[1:]
            
            return jsonify({"smiles": smiles_list})
        else:
            # Default TXT processing - one SMILES per line
            smiles_list = []
            for line in content.splitlines():
                line = line.strip()
                if line:
                    # Take first word as SMILES (handles space-separated data)
                    smiles = line.split()[0]
                    if smiles:
                        smiles_list.append(smiles)
            
            return jsonify({"smiles": smiles_list})
    except Exception as e:
        logger.error(f"Bulk upload error: {str(e)}")
        return jsonify({"error": str(e)}), 500

def _looks_like_smiles(text: str) -> bool:
    """
    Heuristic check to determine if text appears to be a SMILES string.
    
    Uses character analysis and structural indicators to distinguish SMILES
    strings from header text or other data. Not a full validator, but effective
    for filtering out obvious non-SMILES content like headers.
    
    Args:
        text: String to analyze
        
    Returns:
        bool: True if text appears to be SMILES notation
        
    Criteria:
        - Must be a single token (no spaces)
        - >70% of characters must be valid SMILES characters
        - Must contain structural indicators (parentheses, brackets, etc.)
        - Length between 2-999 characters
    """
    # Check for common SMILES characters and reasonable length
    smiles_chars = set("()[]=#CcNnOoSsPpFfClBrI123456789+-@/\\")
    text_chars = set(text)
    
    # If most characters are SMILES-like and it's a single token
    return (
        len(text.split()) == 1 and  # Single token
        len(text_chars.intersection(smiles_chars)) / len(text_chars) > 0.7 and  # Mostly SMILES chars
        any(char in text for char in "()[]=#") and  # Has some structural indicators
        len(text) > 1 and len(text) < 1000  # Reasonable length
    ) 