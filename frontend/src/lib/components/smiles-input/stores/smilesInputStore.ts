import { writable, derived, get } from 'svelte/store';
import { browser } from '$app/environment';
import { appState } from '$lib/stores/appState';

interface SmilesInputState {
  smiles: string;
  isLoading: boolean;
  isFocused: boolean;
  showUploadModal: boolean;
  showGenerateModal: boolean;
  bulkList: string[];
  bulkUploadCallback: ((event: { detail: { list: string[] } }) => void) | null;
}

function createSmilesInputStore() {
  const { subscribe, set, update } = writable<SmilesInputState>({
    smiles: '',
    isLoading: false,
    isFocused: false,
    showUploadModal: false,
    showGenerateModal: false,
    bulkList: [],
    bulkUploadCallback: null,
  });

  return {
    subscribe,
    set,
    update,

    // Smiles text management
    setSmiles: (value: string) => {
      update((s) => ({ ...s, smiles: value }));
      // Sync with global state
      appState.setInputSmiles(value);
    },

    // Load SMILES from global state (for cross-page sync)
    loadFromGlobalState: () => {
      const globalState = get(appState);
      const smilesFromGlobal = globalState.inputSmiles || globalState.currentSmiles;
      if (smilesFromGlobal) {
        update((s) => ({ ...s, smiles: smilesFromGlobal }));
      }
    },

    // Initialize with global state if available
    initWithGlobalState: () => {
      const globalState = get(appState);
      if (globalState.hasFirstPrediction && globalState.currentSmiles) {
        // If there's a previous prediction, load it into input
        update((s) => ({ ...s, smiles: globalState.currentSmiles }));
        appState.setInputSmiles(globalState.currentSmiles);
      }
    },

    // Loading state management
    setLoading: (loading: boolean) => update((s) => ({ ...s, isLoading: loading })),

    // Focus state management
    setFocused: (focused: boolean) => update((s) => ({ ...s, isFocused: focused })),

    // Upload modal management
    showUploadModal: () => {
      if (browser) {
        document.body.classList.add('overlay-open');
      }
      update((s) => ({ ...s, showUploadModal: true }));
    },
    hideUploadModal: () => {
      if (browser) {
        document.body.classList.remove('overlay-open');
      }
      update((s) => ({ ...s, showUploadModal: false }));
    },

    // Generate modal management
    showGenerateModal: () => {
      if (browser) {
        document.body.classList.add('overlay-open');
      }
      update((s) => ({ ...s, showGenerateModal: true }));
    },
    hideGenerateModal: () => {
      if (browser) {
        document.body.classList.remove('overlay-open');
      }
      update((s) => ({ ...s, showGenerateModal: false }));
    },

    // Bulk upload handling
    setBulkUploadCallback: (callback: ((event: { detail: { list: string[] } }) => void) | null) =>
      update((s) => ({ ...s, bulkUploadCallback: callback })),
    triggerBulkUpload: (list: string[]) => {
      let currentState: SmilesInputState | undefined;
      const unsubscribe = subscribe((s) => (currentState = s));
      unsubscribe();

      if (currentState?.bulkUploadCallback) {
        currentState.bulkUploadCallback({ detail: { list } });
      }

      update((s) => ({ ...s, bulkList: list }));
    },

    // Utility methods
    getSmilesList: (): string[] => {
      let currentState: SmilesInputState | undefined;
      const unsubscribe = subscribe((s) => (currentState = s));
      unsubscribe();
      return currentState?.smiles.split(/\s*\n\s*/).filter(Boolean) || [];
    },

    // Clear all data
    reset: () =>
      set({
        smiles: '',
        isLoading: false,
        isFocused: false,
        showUploadModal: false,
        showGenerateModal: false,
        bulkList: [],
        bulkUploadCallback: null,
      }),
  };
}

export const smilesInputStore = createSmilesInputStore();

// Derived stores for computed values
export const lineCount = derived(smilesInputStore, ($store) => $store.smiles.split('\n').length);

export const canSubmit = derived(
  smilesInputStore,
  ($store) => $store.smiles.trim().length > 0 && !$store.isLoading
);
