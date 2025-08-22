"""SMILES generation service."""
from typing import List
from ..client import LLMClient
from ....utils import logger, is_valid_smiles

class SMILESService:
    """
    Service for generating SMILES molecular structures from text descriptions.
    
    Uses LLM to convert natural language descriptions of molecules into
    valid SMILES notation. Includes validation and fallback mechanisms
    to ensure reliable results.
    
    Example descriptions:
    - "aromatic ring with hydroxyl group"
    - "simple alcohol"
    - "caffeine-like structure"
    
    Attributes:
        client: LLMClient for AI-powered generation
    """
    
    def __init__(self):
        self.client = LLMClient()
    
    def generate_smiles(self, description: str, count: int = 1) -> List[str]:
        """
        Generate valid SMILES strings from natural language description.
        
        Converts text descriptions into chemically valid SMILES notation
        using AI, with validation and fallback mechanisms.
        
        Args:
            description: Natural language description of desired molecules
            count: Number of different SMILES to generate
            
        Returns:
            List of valid SMILES strings (up to count requested)
            
        Process:
            1. Generate candidates using google/gemma-3-27b-it:free
            2. Validate each generated SMILES
            3. Use fallbacks if generation fails
            4. Ensure requested count is met when possible
            
        Note:
            Falls back to simple molecules if AI generation fails
        """
        try:
            # Build prompt for SMILES generation
            system_prompt = """You are a chemistry expert. Generate valid SMILES strings only.
Rules:
- Return only SMILES strings, one per line
- No explanations or additional text
- Ensure all SMILES are chemically valid
- Focus on drug-like molecules when possible"""

            user_prompt = f"Generate {count} different SMILES strings for: {description}"
            
            messages = [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ]
            
            response = self.client.chat_completion(
                messages=messages,
                model="google/gemma-3-27b-it:free",
                temperature=0.8  # Higher temperature for more diverse results
            )
            
            # Parse and validate response
            candidates = [line.strip() for line in response.split('\n') if line.strip()]
            valid_smiles = self._validate_smiles_list(candidates)
            
            # Ensure we have enough results
            while len(valid_smiles) < count and len(valid_smiles) < 3:
                # Try again with fallback
                fallback = self._generate_fallback_smiles(count - len(valid_smiles))
                valid_smiles.extend(fallback)
            
            return valid_smiles[:count]
            
        except Exception as e:
            logger.error(f"SMILES generation error: {str(e)}")
            return self._generate_fallback_smiles(count)
    
    def _validate_smiles_list(self, candidates: List[str]) -> List[str]:
        """
        Validate and clean a list of candidate SMILES strings.
        
        Args:
            candidates: Raw strings from LLM response
            
        Returns:
            List of validated SMILES strings
            
        Validation:
            - Extracts first token (handles extra text)
            - Checks chemical validity using RDKit
            - Filters out invalid or empty strings
        """
        valid = []
        for smiles in candidates:
            # Clean up the SMILES string
            cleaned = smiles.strip().split()[0] if smiles.strip() else ""
            if cleaned and is_valid_smiles(cleaned):
                valid.append(cleaned)
        return valid
    
    def _generate_fallback_smiles(self, count: int) -> List[str]:
        """
        Generate fallback SMILES when AI generation fails.
        
        Provides a curated list of simple, common molecules as fallbacks
        to ensure the service always returns valid results.
        
        Args:
            count: Number of SMILES needed
            
        Returns:
            List of simple but valid SMILES strings
            
        Fallback molecules include:
            - Simple alcohols (ethanol, isopropanol)
            - Basic organics (benzene, acetone)
            - Small alkanes (methane, ethane)
        """
        fallbacks = [
            "CCO",  # Ethanol
            "CC(=O)O",  # Acetic acid
            "c1ccccc1",  # Benzene
            "CC(C)O",  # Isopropanol
            "CC(=O)C",  # Acetone
            "CCN",  # Ethylamine
            "CCC",  # Propane
            "CC",   # Ethane
            "C",    # Methane
            "O"     # Water
        ]
        return fallbacks[:count] 