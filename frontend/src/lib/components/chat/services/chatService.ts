import { chatWithSpectrum } from '$lib/services/api';

interface ChatOptions {
  smiles?: string;
  onChunk?: (chunk: string) => void;
}

interface ChatMessage {
  role: 'user' | 'assistant';
  content: string;
}

export const chatService = {
  /**
   * Send a chat message with streaming support
   */
  async sendMessage(messageHistory: ChatMessage[], options: ChatOptions = {}) {
    try {
      const { smiles, onChunk } = options;

      await chatWithSpectrum(messageHistory, {
        smiles,
        onChunk,
      });
    } catch (error) {
      console.error('Chat service error:', error);
      throw new Error(`Failed to send message: ${(error as Error).message}`);
    }
  },

  /**
   * Validate message content
   */
  validateMessage(content: string): boolean {
    if (!content || typeof content !== 'string') return false;
    return content.trim().length > 0;
  },

  /**
   * Sanitize message content
   */
  sanitizeMessage(content: string): string {
    if (!content || typeof content !== 'string') return '';
    return content.trim();
  },

  /**
   * Format error messages for display
   */
  formatErrorMessage(error: unknown): string {
    if (error instanceof Error) {
      return error.message;
    }
    return 'An unexpected error occurred';
  },

  /**
   * Create a chunk processor function
   */
  createChunkProcessor(onMessage: (content: string) => void) {
    let buffer = '';

    return (chunk: string) => {
      if (typeof chunk !== 'string') return;

      buffer += chunk;

      // Process complete sentences for better UX
      const sentences = buffer.split(/[.!?]+/);
      if (sentences.length > 1) {
        const completeSentences = sentences.slice(0, -1).join('. ');
        if (completeSentences.trim()) {
          onMessage(completeSentences + '. ');
          buffer = sentences[sentences.length - 1];
        }
      }

      // Also send chunks periodically for long responses
      if (buffer.length > 100) {
        onMessage(buffer);
        buffer = '';
      }
    };
  },

  /**
   * Validate SMILES format (basic validation)
   */
  isValidSmiles(smiles: string): boolean {
    if (!smiles || typeof smiles !== 'string') return false;

    // Basic SMILES validation - check for valid characters
    const validChars = /^[A-Za-z0-9@+\-\[\]()=#$:.\/\\%]+$/;
    return validChars.test(smiles.trim());
  },

  /**
   * Format timestamp for messages
   */
  formatTimestamp(): string {
    return new Date().toLocaleTimeString();
  },
};
