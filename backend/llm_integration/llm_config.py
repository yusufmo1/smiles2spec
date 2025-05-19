"""
LLM configuration settings for the application.
Defines models, providers, and other LLM-related settings.
"""

import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# OpenRouter base URL
OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1"

# Site information for OpenRouter
SITE_URL = os.environ.get("SITE_URL", "https://smiles2spec.app")
SITE_NAME = os.environ.get("SITE_NAME", "SMILES2Spec App")

# Model configurations
LLM_MODELS = {
    # Primary model
    "primary": "deepseek/deepseek-chat-v3-0324:free",
    
    # Special purpose models
    "chemistry": "deepseek/deepseek-chat-v3-0324:free",  # Best for chemistry tasks
}

# System prompts for different tasks
SYSTEM_PROMPTS = {
    "smiles_generator": "You are a chemical structure generator. Generate valid SMILES strings based on descriptions. Only output the SMILES strings, one per line, without any numbering or prefixes.",
    "general": "You are a helpful assistant specializing in chemistry and molecular structures.",
}

# Helper function to get the appropriate model for a task
def get_model_for_task(task="default"):
    """
    Get the appropriate model for a specific task.
    
    Args:
        task: The task type (default, chemistry)
        
    Returns:
        Model ID string
    """
    if task == "chemistry" and "chemistry" in LLM_MODELS:
        return LLM_MODELS["chemistry"]
    else:
        return LLM_MODELS["primary"]

# Helper function to get the system prompt for a task
def get_system_prompt(task_type):
    """
    Get the system prompt for a specific task type.
    
    Args:
        task_type: Type of task (smiles_generator, general, etc.)
        
    Returns:
        System prompt string
    """
    return SYSTEM_PROMPTS.get(task_type, SYSTEM_PROMPTS["general"])

# OpenRouter API request helper
def get_openrouter_headers():
    """
    Get the standard headers for OpenRouter API requests.
    
    Returns:
        Dictionary of headers
    """
    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        raise ValueError("OPENROUTER_API_KEY environment variable is not set")
        
    return {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}",
        "HTTP-Referer": SITE_URL,
        "X-Title": SITE_NAME,
    } 