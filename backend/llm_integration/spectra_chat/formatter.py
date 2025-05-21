"""Format spectrum data for inclusion in chat."""
import base64
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np


def spectrum_to_png_base64(spectrum_data):
    """Convert spectrum data to a PNG image encoded as base64."""
    # Extract spectrum data
    if not spectrum_data or 'spectrum' not in spectrum_data:
        return None

    try:
        x = spectrum_data['spectrum']['x']
        y = spectrum_data['spectrum']['y']

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
        return png_base64

    except Exception as e:
        print(f"Error creating spectrum image: {e}")
        return None


def create_spectrum_message(spectrum_data):
    """Create a message object with the spectrum image."""
    png_base64 = spectrum_to_png_base64(spectrum_data)

    if not png_base64:
        return None

    # Create message with image
    return {
        "role": "user",
        "content": [
            {
                "type": "text",
                "text": "Here is the predicted mass spectrum for the molecule:",
            },
            {
                "type": "image_url",
                "image_url": {
                    "url": f"data:image/png;base64,{png_base64}"
                }
            }
        ]
    }
