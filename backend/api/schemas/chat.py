"""Chat API schemas."""
from pydantic import BaseModel, validator
from typing import List, Optional, Dict, Any

class ChatMessage(BaseModel):
    """Chat message schema."""
    role: str
    content: str
    
    @validator('role')
    def validate_role(cls, v):
        allowed_roles = {'user', 'assistant', 'system'}
        if v not in allowed_roles:
            raise ValueError(f'Role must be one of: {allowed_roles}')
        return v
    
    @validator('content')
    def validate_content(cls, v):
        # Convert to string and strip
        content = str(v).strip() if v is not None else ""
        
        # Allow empty content for system messages or in certain contexts
        if not content:
            return ""  # Return empty string instead of raising error
            
        if len(content) > 10000:
            raise ValueError('Message content too long (max 10000 characters)')
        return content

class ChatRequest(BaseModel):
    """Chat request schema."""
    messages: List[ChatMessage]
    stream: Optional[bool] = False
    smiles: Optional[str] = None
    
    @validator('messages')
    def validate_messages(cls, v):
        if not v:
            raise ValueError('Messages list cannot be empty')
        if len(v) > 50:
            raise ValueError('Maximum 50 messages allowed')
            
        # Filter out empty messages
        valid_messages = [msg for msg in v if msg.content.strip()]
        if not valid_messages:
            raise ValueError('At least one message must have content')
            
        return valid_messages  # Return filtered messages
    
    @validator('smiles')
    def validate_smiles(cls, v):
        if v is not None and (not str(v).strip() or len(str(v).strip()) > 1000):
            raise ValueError('Invalid SMILES string')
        return str(v).strip() if v else None

class ChatResponse(BaseModel):
    """Chat response schema."""
    message: str
    metadata: Optional[Dict[str, Any]] = None

class GenerateSMILESRequest(BaseModel):
    """Generate SMILES request schema."""
    count: Optional[int] = 1
    description: Optional[str] = ""
    
    @validator('count')
    def validate_count(cls, v):
        if v < 1 or v > 10:
            raise ValueError('Count must be between 1 and 10')
        return v
    
    @validator('description')
    def validate_description(cls, v):
        if v and len(v) > 500:
            raise ValueError('Description too long (max 500 characters)')
        return v or ""

class GenerateSMILESResponse(BaseModel):
    """Generate SMILES response schema."""
    smiles: List[str]
    metadata: Optional[Dict[str, Any]] = None 