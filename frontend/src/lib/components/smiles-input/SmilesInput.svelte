<script lang="ts">
  /**
   * SmilesInput Component
   *
   * Main input interface for entering SMILES molecular structures.
   * Provides multiple input methods including direct typing, file upload,
   * and AI-powered generation from descriptions.
   *
   * Features:
   * - Multi-line SMILES input with auto-detection of bulk mode
   * - File upload support (CSV, TXT)
   * - AI generation from natural language
   * - Visual feedback for focus state
   * - Responsive design for mobile
   * - Cross-page state synchronization
   *
   * Events:
   * - predict: Single SMILES prediction
   * - bulk: Multiple SMILES batch processing
   *
   * @component
   */
  import { createEventDispatcher, onMount } from 'svelte';
  import { smilesInputStore } from './stores/smilesInputStore';
  import SmilesTextArea from './components/SmilesTextArea.svelte';
  import ActionButtons from './components/ActionButtons.svelte';
  import UploadModal from './components/UploadModal.svelte';
  import GenerateModal from './components/GenerateModal.svelte';

  const dispatch = createEventDispatcher();

  /**
   * Handle form submission
   * Determines whether to process as single or bulk based on line count
   */
  function handleSubmit() {
    if (!$smilesInputStore.smiles.trim() || $smilesInputStore.isLoading) return;

    smilesInputStore.setLoading(true);

    // Get all non-empty lines
    const lines = smilesInputStore.getSmilesList();

    // If there's only one line, treat as a single SMILES
    if (lines.length === 1) {
      dispatch('predict', { smiles: lines[0] });
    } else if (lines.length > 1) {
      // If there are multiple lines, send as bulk
      dispatch('bulk', { list: lines });
    }
  }

  /**
   * Show file upload modal
   */
  function handleUpload() {
    smilesInputStore.showUploadModal();
  }

  /**
   * Show AI generation modal
   */
  function handleGenerate() {
    smilesInputStore.showGenerateModal();
  }

  /**
   * Set loading state - exposed for parent components
   * @public
   */
  export function setLoading(loading: boolean) {
    smilesInputStore.setLoading(loading);
  }

  // Initialize with global state on mount
  onMount(() => {
    smilesInputStore.initWithGlobalState();
  });
</script>

<div class="input-container glass-card" class:focused={$smilesInputStore.isFocused}>
  <SmilesTextArea on:submit={handleSubmit} />

  <ActionButtons on:submit={handleSubmit} on:upload={handleUpload} on:generate={handleGenerate} />
</div>

<!-- Modals -->
<UploadModal />
<GenerateModal />

<style>
  .input-container {
    display: flex;
    align-items: center;
    margin: 0 0 2rem;
    padding: 0.65rem;
    transition: all var(--transition-elastic);
    background: var(--surface-glass-hover);
    border-radius: var(--enforce-pill);
    position: relative;
    overflow: hidden !important;
  }

  .input-container.focused {
    box-shadow:
      0 10px 35px 0 rgba(0, 0, 0, 0.08),
      0 0 0 2px var(--accent-soft),
      0 2px 8px 0 rgba(255, 255, 255, 0.4) inset;
  }

  /* Highlight effect for focused state */
  .input-container.focused::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 40%;
    background: linear-gradient(to bottom, rgba(255, 255, 255, 0.25), rgba(255, 255, 255, 0.05));
    z-index: -1;
  }

  @media (max-width: 640px) {
    .input-container {
      flex-direction: column;
      padding: 1.25rem;
    }
  }
</style>
