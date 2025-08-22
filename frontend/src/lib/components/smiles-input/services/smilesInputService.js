import { uploadSmilesFile, generateRandomSmiles } from '$lib/services/api';

export const smilesInputService = {
  /**
   * Process uploaded file and extract SMILES data
   * @param {File} file - The uploaded file
   * @returns {Promise<{smiles: string[]}>} - Processed SMILES data
   */
  async processFile(file) {
    const validTypes = ['text/plain', 'text/csv', 'application/csv', 'application/vnd.ms-excel'];
    const fileExt = file.name.split('.').pop().toLowerCase();

    if (!file || !(validTypes.includes(file.type) || ['txt', 'csv'].includes(fileExt))) {
      throw new Error('Please upload a text (.txt) or CSV (.csv) file');
    }

    try {
      const result = await uploadSmilesFile(file);
      return result;
    } catch (error) {
      throw new Error(`Failed to process file: ${error.message}`);
    }
  },

  /**
   * Generate random SMILES strings
   * @param {Object} options - Generation options
   * @param {number} options.count - Number of SMILES to generate
   * @param {string} options.description - Description for AI generation
   * @returns {Promise<{smiles: string[]}>} - Generated SMILES data
   */
  async generateSmiles(options) {
    try {
      const result = await generateRandomSmiles(options.count, options.description);
      return { smiles: result.smiles };
    } catch (error) {
      throw new Error(`Failed to generate SMILES: ${error.message}`);
    }
  },

  /**
   * Validate SMILES string format
   * @param {string} smiles - SMILES string to validate
   * @returns {boolean} - Whether the SMILES is valid
   */
  validateSmiles(smiles) {
    if (!smiles || typeof smiles !== 'string') {
      return false;
    }

    const trimmed = smiles.trim();
    if (trimmed.length === 0) {
      return false;
    }

    // Basic SMILES validation - check for common invalid characters
    const invalidChars = /[^A-Za-z0-9\[\]()=#@+\-\.\\\/]/;
    return !invalidChars.test(trimmed);
  },

  /**
   * Parse multiple SMILES from text input
   * @param {string} text - Multi-line SMILES text
   * @returns {string[]} - Array of individual SMILES strings
   */
  parseSmilesList(text) {
    if (!text || typeof text !== 'string') {
      return [];
    }

    return text
      .split(/\s*\n\s*/)
      .map((line) => line.trim())
      .filter((line) => line.length > 0 && this.validateSmiles(line));
  },

  /**
   * Format SMILES list for display
   * @param {string[]} smilesList - Array of SMILES strings
   * @returns {string} - Formatted multi-line string
   */
  formatSmilesForDisplay(smilesList) {
    if (!Array.isArray(smilesList)) {
      return '';
    }

    return smilesList.filter((smiles) => this.validateSmiles(smiles)).join('\n');
  },
};
