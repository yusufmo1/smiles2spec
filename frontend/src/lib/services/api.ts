import { browser } from '$app/environment';

/**
 * API Service Module
 *
 * Provides a comprehensive TypeScript client for interacting with the
 * SMILES2SPEC backend API. Handles all HTTP communication, error handling,
 * and data transformation between the frontend and backend services.
 *
 * Features:
 * - Automatic environment detection (dev/prod)
 * - Type-safe interfaces for all API responses
 * - Streaming support for chat functionality
 * - Error handling with descriptive messages
 * - Browser-only execution guards
 */

// API Configuration - Updated for subdomain in production
function getApiUrl(): string {
  if (!browser) return 'http://localhost:5050';

  const hostname = window.location.hostname;

  // Development environment - use proxy
  if (hostname === 'localhost' || hostname === '127.0.0.1') {
    return '/api'; // Proxied to localhost:5050 via Vite
  }

  // Production environment - use subdomain
  return 'https://api.spectralsimulation.com';
}

// TypeScript interfaces

/**
 * Spectrum data structure containing x (m/z) and y (intensity) arrays
 */
export interface SpectrumData {
  x: number[];
  y: number[];
}

/**
 * Individual peak in a mass spectrum
 */
export interface PeakData {
  mz: number;
  intensity: number;
}

/**
 * Complete response from spectrum prediction endpoint
 */
export interface PredictionResponse {
  spectrum: SpectrumData;
  peaks: PeakData[];
  structure_png: string;
  chemical_name: string;
  molecular_weight?: number;
  exact_mass?: number;
  smiles: string;
}

/**
 * Response from SMILES generation endpoint
 */
export interface GeneratedSmilesResponse {
  smiles: string[];
  success: boolean;
  message?: string;
}

/**
 * Response from bulk SMILES file upload
 */
export interface BulkUploadResponse {
  smiles: string[];
  success: boolean;
  message?: string;
}

/**
 * Chat message format for AI conversations
 */
export interface ChatMessage {
  role: string;
  content: string;
}

/**
 * Options for chat API calls
 */
export interface ChatOptions {
  smiles?: string;
  onChunk?: (chunk: string) => void;
}

/**
 * Predict mass spectrum from SMILES molecular structure
 *
 * Sends a SMILES string to the backend ML model and receives predicted
 * mass spectrum data including peaks, molecular properties, and structure
 * visualization.
 *
 * @param smiles - Valid SMILES string representation of molecule
 * @returns Promise resolving to complete prediction data
 * @throws Error if SMILES is invalid, server unreachable, or prediction fails
 *
 * @example
 * const result = await predictSpectrum('CCO'); // Predict ethanol spectrum
 * console.log(result.chemical_name); // "Ethanol"
 * console.log(result.spectrum.x); // m/z values
 * console.log(result.spectrum.y); // intensities
 */
export async function predictSpectrum(smiles: string): Promise<PredictionResponse> {
  if (!browser) {
    throw new Error('API calls only available in browser');
  }

  if (!smiles.trim()) {
    throw new Error('SMILES string is required');
  }

  const API_URL = getApiUrl();

  try {
    const response = await fetch(`${API_URL}/predict`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        Accept: 'application/json',
      },
      body: JSON.stringify({ smiles: smiles.trim() }),
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Prediction failed (${response.status}): ${errorText}`);
    }

    const data = await response.json();

    // Validate response structure
    if (!data.spectrum || !data.spectrum.x || !data.spectrum.y) {
      throw new Error('Invalid response: missing spectrum data');
    }

    return {
      spectrum: data.spectrum,
      peaks: data.peaks || [],
      structure_png: data.structure_png || '',
      chemical_name: data.chemical_name || '',
      molecular_weight: data.molecular_weight,
      exact_mass: data.exact_mass,
      smiles: smiles,
    };
  } catch (error) {
    console.error('Prediction error:', error);

    if (error instanceof TypeError && error.message.includes('fetch')) {
      throw new Error('Cannot connect to prediction server. Please ensure the backend is running.');
    }

    throw error;
  }
}

/**
 * Generate SMILES strings from natural language descriptions
 *
 * Uses AI to generate chemically valid SMILES strings based on text
 * descriptions. Useful for exploring molecular structures without
 * knowing exact SMILES notation.
 *
 * @param count - Number of SMILES to generate (1-10)
 * @param description - Optional text description of desired molecules
 * @returns Promise resolving to array of generated SMILES
 * @throws Error if generation fails or server unreachable
 *
 * @example
 * // Generate 3 aromatic compounds
 * const result = await generateRandomSmiles(3, 'aromatic ring compounds');
 * console.log(result.smiles); // ['c1ccccc1', 'c1ccc(O)cc1', ...]
 */
export async function generateRandomSmiles(
  count: number = 1,
  description?: string
): Promise<GeneratedSmilesResponse> {
  if (!browser) {
    throw new Error('API calls only available in browser');
  }

  const API_URL = getApiUrl();

  try {
    const response = await fetch(`${API_URL}/generate_smiles`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        Accept: 'application/json',
      },
      body: JSON.stringify({
        count: Math.max(1, Math.min(count, 10)), // Limit between 1-10
        description: description?.trim() || undefined,
      }),
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Generation failed (${response.status}): ${errorText}`);
    }

    return response.json();
  } catch (error) {
    console.error('SMILES generation error:', error);

    if (error instanceof TypeError && error.message.includes('fetch')) {
      throw new Error('Cannot connect to generation server. Please ensure the backend is running.');
    }

    throw error;
  }
}

/**
 * Export spectrum data as MSP format file
 *
 * Generates MSP (Mass Spectral Peak) format file for a single molecule.
 * MSP is a standard format for mass spectrometry data exchange.
 *
 * @param smiles - SMILES string to export
 * @returns Promise resolving to MSP file as Blob
 * @throws Error if export fails
 *
 * @example
 * const blob = await exportMsp('CCO');
 * const url = URL.createObjectURL(blob);
 * // Trigger download...
 */
export async function exportMsp(smiles: string): Promise<Blob> {
  if (!browser) {
    throw new Error('API calls only available in browser');
  }

  const API_URL = getApiUrl();

  try {
    const response = await fetch(`${API_URL}/export_msp`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles }),
    });

    if (!response.ok) {
      throw new Error(`MSP export failed (${response.status})`);
    }

    return response.blob();
  } catch (error) {
    console.error('MSP export error:', error);
    throw error;
  }
}

/**
 * Export multiple SMILES as batch MSP file
 *
 * Generates a single MSP file containing spectra for multiple molecules.
 * Useful for creating spectral libraries or bulk analysis.
 *
 * @param smilesList - Array of SMILES strings to export
 * @returns Promise resolving to combined MSP file as Blob
 * @throws Error if batch export fails
 *
 * @example
 * const smiles = ['CCO', 'CC(=O)O', 'c1ccccc1'];
 * const blob = await exportMspBatch(smiles);
 */
export async function exportMspBatch(smilesList: string[]): Promise<Blob> {
  if (!browser) {
    throw new Error('API calls only available in browser');
  }

  const API_URL = getApiUrl();

  try {
    const response = await fetch(`${API_URL}/export_msp_batch`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles_list: smilesList }),
    });

    if (!response.ok) {
      throw new Error(`Batch MSP export failed (${response.status})`);
    }

    return response.blob();
  } catch (error) {
    console.error('Batch MSP export error:', error);
    throw error;
  }
}

/**
 * Upload SMILES file for bulk processing
 *
 * Processes CSV or text files containing multiple SMILES strings.
 * Supports both comma-separated and line-separated formats.
 *
 * @param file - File object containing SMILES data
 * @returns Promise resolving to extracted SMILES array
 * @throws Error if file processing fails
 *
 * @example
 * const file = new File(['CCO\nCC(=O)O'], 'molecules.txt');
 * const result = await uploadSmilesFile(file);
 * console.log(result.smiles); // ['CCO', 'CC(=O)O']
 */
export async function uploadSmilesFile(file: File): Promise<BulkUploadResponse> {
  if (!browser) {
    throw new Error('API calls only available in browser');
  }

  const API_URL = getApiUrl();

  try {
    const formData = new FormData();
    formData.append('file', file);

    const response = await fetch(`${API_URL}/smiles_bulk`, {
      method: 'POST',
      body: formData,
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`File upload failed (${response.status}): ${errorText}`);
    }

    return response.json();
  } catch (error) {
    console.error('File upload error:', error);

    if (error instanceof TypeError && error.message.includes('fetch')) {
      throw new Error('Cannot connect to upload server. Please ensure the backend is running.');
    }

    throw error;
  }
}

/**
 * Chat with AI about mass spectrometry using streaming responses
 *
 * Enables real-time conversation with an AI assistant specialized in
 * mass spectrometry. Supports context-aware discussions when a SMILES
 * string is provided. Uses server-sent events for streaming responses.
 *
 * @param messages - Chat history with role/content pairs
 * @param options - Optional configuration including SMILES context and chunk handler
 * @returns Promise that resolves when streaming completes
 * @throws Error if chat fails or connection lost
 *
 * @example
 * await chatWithSpectrum(
 *   [{role: 'user', content: 'Explain this spectrum'}],
 *   {
 *     smiles: 'CCO',
 *     onChunk: (chunk) => console.log(chunk)
 *   }
 * );
 */
export async function chatWithSpectrum(
  messages: ChatMessage[],
  options: ChatOptions = {}
): Promise<void> {
  if (!browser) {
    throw new Error('API calls only available in browser');
  }

  const API_URL = getApiUrl();
  const { smiles, onChunk } = options;

  try {
    const payload: any = {
      messages,
      stream: true,
    };

    if (smiles) {
      payload.smiles = smiles;
    }

    const response = await fetch(`${API_URL}/chat`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json; charset=utf-8',
        Accept: 'text/event-stream; charset=utf-8',
      },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Chat failed (${response.status}): ${errorText}`);
    }

    // Handle streaming response
    const reader = response.body!.pipeThrough(new TextDecoderStream('utf-8')).getReader();

    let buffer = '';

    try {
      while (true) {
        const { done, value } = await reader.read();

        if (done) {
          // Process any remaining buffer content
          if (buffer.trim() && typeof onChunk === 'function') {
            try {
              const parsed = JSON.parse(buffer.trim());
              if (parsed.chunk) {
                onChunk(parsed.chunk);
              }
            } catch (e) {
              console.warn('Error parsing final buffer:', e);
            }
          }
          break;
        }

        buffer += value;

        // Process complete lines from buffer
        let lineEnd;
        while ((lineEnd = buffer.indexOf('\n')) !== -1) {
          const line = buffer.slice(0, lineEnd).trim();
          buffer = buffer.slice(lineEnd + 1);

          if (line.startsWith('data: ')) {
            const data = line.slice(6);

            if (data === '[DONE]') {
              return;
            }

            try {
              const parsed = JSON.parse(data);
              if (parsed.chunk && typeof onChunk === 'function') {
                onChunk(parsed.chunk);
              }
            } catch (e) {
              console.warn('Error parsing chunk:', e, 'Raw data:', data);
            }
          }
        }
      }
    } finally {
      reader.releaseLock();
    }
  } catch (error) {
    console.error('Chat error:', error);

    if (error instanceof TypeError && error.message.includes('fetch')) {
      throw new Error('Cannot connect to chat server. Please ensure the backend is running.');
    }

    throw error;
  }
}

/**
 * Health check for the API server
 *
 * Verifies that the backend API is running and responsive.
 * Used for connection status indicators and error handling.
 *
 * @returns Promise resolving to true if API is healthy, false otherwise
 *
 * @example
 * const isHealthy = await checkApiHealth();
 * if (!isHealthy) {
 *   console.error('Backend API is not responding');
 * }
 */
export async function checkApiHealth(): Promise<boolean> {
  if (!browser) return false;

  const API_URL = getApiUrl();

  try {
    const response = await fetch(`${API_URL}/health`, {
      method: 'GET',
      headers: { Accept: 'application/json' },
    });

    return response.ok;
  } catch (error) {
    console.warn('API health check failed:', error);
    return false;
  }
}
