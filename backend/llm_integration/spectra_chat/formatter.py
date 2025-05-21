"""Format spectrum data for inclusion in chat."""
import base64
from io import BytesIO
import logging
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

def spectrum_to_png_base64(spectrum_data):
    """Convert spectrum data to a PNG image encoded as base64."""
    # Extract spectrum data
    if not spectrum_data or 'spectrum' not in spectrum_data:
        logger.warning("No spectrum data available for visualization")
        return None
        
    try:
        x = spectrum_data['spectrum'].get('x', [])
        y = spectrum_data['spectrum'].get('y', [])
        
        if not x or not y or len(x) != len(y):
            logger.warning(f"Invalid spectrum data: x={len(x)} points, y={len(y)} points")
            return None
        
        # Create plot
        plt.figure(figsize=(10, 6))
        plt.bar(x, y, width=0.5, alpha=0.7)
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title(f"Predicted Spectrum: {spectrum_data.get('chemical_name', spectrum_data.get('smiles', 'Unknown'))}")
        
        # Convert plot to PNG
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=100)
        plt.close()
        buf.seek(0)
        
        # Encode as base64
        png_base64 = base64.b64encode(buf.read()).decode('utf-8')
        logger.info(f"Successfully generated spectrum image ({len(png_base64)} bytes)")
        return png_base64
        
    except Exception as e:
        logger.error(f"Error creating spectrum image: {e}", exc_info=True)
        return None

def create_spectrum_message(spectrum_data):
    """Create a message object with the spectrum image."""
    try:
        png_base64 = spectrum_to_png_base64(spectrum_data)
        
        if not png_base64:
            logger.warning("Could not create spectrum image, skipping visual context")
            return None
            
        # Create message with image
        return {
            "role": "user",
            "content": [
                {
                    "type": "text",
                    "text": "Here is the predicted mass spectrum for the molecule:"
                },
                {
                    "type": "image_url",
                    "image_url": {
                        "url": f"data:image/png;base64,{png_base64}"
                    }
                }
            ]
        }
    except Exception as e:
        logger.error(f"Error creating spectrum message: {e}", exc_info=True)
        return None
