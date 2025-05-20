"""
Module for downloading ML models from Azure Blob Storage.
"""

import os
from pathlib import Path
from azure.storage.blob import BlobServiceClient
from azure.core.exceptions import ResourceNotFound

def ensure_model_directory(model_path: str) -> None:
    """Ensure the model directory exists."""
    Path(os.path.dirname(model_path)).mkdir(parents=True, exist_ok=True)

def download_spectrum_predictor() -> None:
    """
    Download the spectrum predictor model from Azure Blob Storage if it doesn't exist locally.
    Uses environment variables for configuration.
    """
    connection_string = os.getenv('AZURE_STORAGE_CONNECTION_STRING')
    container_name = os.getenv('AZURE_STORAGE_CONTAINER')
    blob_name = os.getenv('AZURE_STORAGE_BLOB')
    local_path = os.getenv('LOCAL_MODEL_PATH')

    if not all([connection_string, container_name, blob_name, local_path]):
        raise ValueError("Missing required environment variables for model download")

    # If model already exists locally, skip download
    if os.path.exists(local_path):
        print(f"Model already exists at {local_path}")
        return

    print(f"Downloading model from Azure Blob Storage...")
    
    try:
        # Create the BlobServiceClient
        blob_service_client = BlobServiceClient.from_connection_string(connection_string)
        
        # Get container client
        container_client = blob_service_client.get_container_client(container_name)
        
        # Get blob client
        blob_client = container_client.get_blob_client(blob_name)

        # Ensure the model directory exists
        ensure_model_directory(local_path)

        # Download the blob
        with open(local_path, "wb") as model_file:
            download_stream = blob_client.download_blob()
            model_file.write(download_stream.readall())

        print(f"Successfully downloaded model to {local_path}")

    except ResourceNotFound:
        raise Exception(f"Model file {blob_name} not found in Azure Storage container {container_name}")
    except Exception as e:
        raise Exception(f"Failed to download model: {str(e)}")

def initialize_models():
    """
    Initialize all required models by downloading them if necessary.
    Currently only downloads the spectrum predictor model.
    """
    download_spectrum_predictor() 