/**
 * Global Application State Store
 *
 * Centralized state management for the SMILES2SPEC application using
 * Svelte stores. Maintains all application-wide data including current
 * predictions, bulk processing state, and console logging.
 *
 * This store serves as the single source of truth for:
 * - Current molecular prediction data
 * - Bulk SMILES processing queue
 * - Application loading/error states
 * - Console output for user feedback
 */

import { writable, derived, get } from 'svelte/store';
import type { PredictionResponse } from '$lib/services/api';

/**
 * Main application state interface
 */
export interface GlobalAppState {
  // Current prediction data
  currentSmiles: string;
  currentName: string;
  spectrumData: any | null;
  peakData: any[];
  structurePNG: string | null;
  molecularWeight: number | null;
  exactMass: number | null;

  // Input state - tracks what user is currently typing
  inputSmiles: string;

  // Bulk processing
  bulkList: string[];
  bulkIndex: number;

  // Prediction state
  hasFirstPrediction: boolean;
  isLoading: boolean;
  error: string | null;

  // Console data
  consoleEntries: ConsoleEntry[];
}

/**
 * Console log entry for activity tracking
 */
export interface ConsoleEntry {
  timestamp: Date;
  type: 'info' | 'success' | 'error' | 'cpu';
  message: string;
  details?: Record<string, any>;
}

const initialState: GlobalAppState = {
  currentSmiles: '',
  currentName: '',
  spectrumData: null,
  peakData: [],
  structurePNG: null,
  molecularWeight: null,
  exactMass: null,
  inputSmiles: '',
  bulkList: [],
  bulkIndex: 0,
  hasFirstPrediction: false,
  isLoading: false,
  error: null,
  consoleEntries: [
    {
      timestamp: new Date(),
      type: 'info',
      message: 'System initialized. Ready for spectrum prediction.',
      details: {},
    },
  ],
};

/**
 * Create the main application state store with methods
 *
 * Provides a comprehensive API for managing application state including
 * prediction data, bulk processing, console logging, and navigation.
 *
 * @returns Svelte store with custom methods for state management
 */
function createAppState() {
  const { subscribe, set, update } = writable<GlobalAppState>(initialState);

  return {
    subscribe,

    /**
     * Set prediction data from API response
     * Updates all relevant fields and marks first prediction complete
     */
    setPredictionData: (data: PredictionResponse) => {
      update((state) => ({
        ...state,
        currentSmiles: data.smiles,
        currentName: data.chemical_name || '',
        spectrumData: data.spectrum,
        peakData: data.peaks || [],
        structurePNG: data.structure_png || null,
        molecularWeight: data.molecular_weight || null,
        exactMass: data.exact_mass || null,
        inputSmiles: data.smiles, // Sync input with prediction
        hasFirstPrediction: true,
        error: null,
      }));
    },

    /**
     * Set the current input SMILES (what user is typing)
     * This allows syncing across pages
     */
    setInputSmiles: (smiles: string) => {
      update((state) => ({ ...state, inputSmiles: smiles }));
    },

    /**
     * Load a previous prediction into the input
     * Useful for cross-page navigation
     */
    loadPredictionToInput: () => {
      update((state) => ({
        ...state,
        inputSmiles: state.currentSmiles,
      }));
    },

    /**
     * Set list of SMILES for bulk processing
     * Resets index to beginning of list
     */
    setBulkList: (list: string[]) => {
      update((state) => ({
        ...state,
        bulkList: list,
        bulkIndex: 0,
      }));
    },

    /**
     * Navigate to specific index in bulk list
     * Clamps to valid range [0, length-1]
     */
    setBulkIndex: (index: number) => {
      update((state) => ({
        ...state,
        bulkIndex: Math.max(0, Math.min(index, state.bulkList.length - 1)),
      }));
    },

    /**
     * Add new entry to console log
     * Automatically timestamps the entry
     */
    addConsoleEntry: (entry: Omit<ConsoleEntry, 'timestamp'>) => {
      update((state) => ({
        ...state,
        consoleEntries: [
          ...state.consoleEntries,
          {
            ...entry,
            timestamp: new Date(),
          },
        ],
      }));
    },

    /**
     * Set global loading state
     */
    setLoading: (loading: boolean) => {
      update((state) => ({ ...state, isLoading: loading }));
    },

    /**
     * Set global error message
     */
    setError: (error: string | null) => {
      update((state) => ({ ...state, error }));
    },

    /**
     * Get current SMILES based on bulk processing state
     * Returns bulk item if processing, otherwise current SMILES
     */
    getCurrentSmiles: () => {
      const state = get({ subscribe });
      return state.bulkList.length > 0 ? state.bulkList[state.bulkIndex] : state.currentSmiles;
    },

    /**
     * Reset entire application state to initial values
     */
    reset: () => set(initialState),
  };
}

/**
 * Main application state store instance
 */
export const appState = createAppState();

/**
 * Derived store for current prediction data
 * Provides convenient access to all prediction-related fields
 */
export const currentPredictionData = derived(appState, ($state) => ({
  smiles: $state.currentSmiles,
  name: $state.currentName,
  spectrum: $state.spectrumData,
  peaks: $state.peakData,
  structure: $state.structurePNG,
  molecularWeight: $state.molecularWeight,
  exactMass: $state.exactMass,
}));

/**
 * Derived store for formatted console text
 * Converts console entries to timestamped text format
 */
export const consoleText = derived(appState, ($state) =>
  $state.consoleEntries
    .map((entry) => `[${entry.timestamp.toLocaleTimeString()}] ${entry.message}`)
    .join('\n')
);
