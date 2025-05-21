#!/usr/bin/env python
"""
Test script for SMILES generator.
"""
import os
import logging
import sys
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(level=logging.DEBUG, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def test_smiles_generator():
    """Test the SMILES generator with a simple description."""
    try:
        logger.info("Testing SMILES generator...")
        
        # Print environment variables for debugging
        api_key = os.environ.get("OPENROUTER_API_KEY")
        logger.info(f"API key present: {bool(api_key)}")
        if api_key:
            logger.info(f"API key prefix: {api_key[:5]}...")
        
        # Import the module
        logger.info("Importing module...")
        from llm_integration.smiles_generator.service import generate_smiles
        
        # Test the function
        logger.info("Calling generate_smiles...")
        result = generate_smiles(
            description="molecule with a benzene ring", 
            count=2
        )
        logger.info(f"Result: {result}")
        return result
    except Exception as e:
        logger.error(f"Error: {str(e)}", exc_info=True)
        return None

def test_module_structure():
    """Test the module structure to ensure imports are working."""
    try:
        # Root module
        logger.info("Importing llm_integration...")
        import llm_integration
        logger.info(f"llm_integration.__all__: {llm_integration.__all__}")
        
        # Import the generate_smiles function from root
        logger.info("Importing generate_smiles from root...")
        from llm_integration import generate_smiles
        logger.info(f"generate_smiles function: {generate_smiles}")
        
        # Test the function from root import
        logger.info("Calling generate_smiles from root import...")
        result = generate_smiles(
            description="simple alcohol", 
            count=1
        )
        logger.info(f"Result from root import: {result}")
        return result
    except Exception as e:
        logger.error(f"Error in module structure test: {str(e)}", exc_info=True)
        return None

if __name__ == "__main__":
    # First test the module structure
    logger.info("=== Testing module structure ===")
    test_module_structure()
    
    # Then test the direct import
    logger.info("\n=== Testing direct import ===")
    test_smiles_generator() 