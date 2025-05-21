"""Shared configuration for LLM integration."""
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# API configurations
OPENROUTER_API_KEY = os.environ.get("OPENROUTER_API_KEY")
OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1"

# Site information
SITE_URL = os.environ.get("SITE_URL", "https://smiles2spec.app")
SITE_NAME = os.environ.get("SITE_NAME", "SMILES2Spec App")

# Model configurations - using free tier models available on OpenRouter
MODELS = {
    "chat": "google/gemma-3-27b-it:free",    # General chat, free tier
    "chemistry": "google/gemma-3-27b-it:free",  # For chemistry tasks
}

def get_headers():
    """Get standard headers for API requests."""
    if not OPENROUTER_API_KEY:
        raise ValueError("OPENROUTER_API_KEY is not set")
        
    return {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {OPENROUTER_API_KEY}",
        "HTTP-Referer": SITE_URL,
        "X-Title": SITE_NAME,
    }
