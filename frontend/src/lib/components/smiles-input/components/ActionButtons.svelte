<script>
  import { createEventDispatcher } from 'svelte';
  import { smilesInputStore, canSubmit } from '../stores/smilesInputStore';
  import { DynamicIcon } from '$lib/components/icons';

  const dispatch = createEventDispatcher();

  function handleSubmit() {
    dispatch('submit');
  }

  function handleUpload() {
    dispatch('upload');
  }

  function handleGenerate() {
    dispatch('generate');
  }
</script>

<div class="action-buttons-container">
  <!-- Magic Wand button -->
  <button class="pill-button generate" aria-label="Generate SMILES" on:click={handleGenerate}>
    <DynamicIcon name="WandIcon" size={18} color="white" />
  </button>

  <!-- Upload button -->
  <button class="pill-button upload" aria-label="Upload SMILES list" on:click={handleUpload}>
    <DynamicIcon name="UploadIcon" size={18} color="white" />
  </button>

  <!-- Analyze button -->
  <button
    class="pill-button"
    on:click={handleSubmit}
    disabled={!$canSubmit}
    aria-label="Analyze molecule"
  >
    {#if $smilesInputStore.isLoading}
      <DynamicIcon name="LoadingIcon" size={16} />
      <span>Processing</span>
    {:else}
      <span>Analyze</span>
    {/if}
  </button>
</div>

<style>
  .action-buttons-container {
    display: flex;
    align-items: center;
    gap: 0.25rem;
  }

  .pill-button {
    padding: 0.85rem 1.8rem;
    border-radius: var(--enforce-pill);
    transition: all var(--transition-elastic);
    font-weight: 500;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    color: white;
    border: none;
    font-size: 0.95rem;
    cursor: pointer;
    box-shadow:
      0 6px 18px rgba(120, 121, 255, 0.25),
      0 0 0 1px rgba(255, 255, 255, 0.1);
  }

  .generate {
    background: #9c27b0; /* Purple color for wand button */
    padding: 0.65rem;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .generate:hover {
    background: #7b1fa2;
  }

  .upload {
    background: var(--accent-secondary);
    padding: 0.65rem;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .upload:hover {
    background: var(--accent);
  }

  .pill-button:hover:not(:disabled) {
    background: linear-gradient(135deg, var(--accent-secondary) 0%, var(--accent) 100%);
    transform: translateY(-2px) scale(1.025);
    box-shadow:
      0 8px 24px rgba(120, 121, 255, 0.35),
      0 0 0 1px rgba(255, 255, 255, 0.15);
  }

  .pill-button:disabled {
    background: rgba(0, 0, 0, 0.1);
    color: var(--text-tertiary);
    cursor: not-allowed;
    box-shadow: none;
  }

  @keyframes spin {
    0% {
      transform: rotate(0deg);
    }
    100% {
      transform: rotate(360deg);
    }
  }

  @media (max-width: 640px) {
    .action-buttons-container {
      width: 100%;
      justify-content: space-between;
      gap: 0.5rem;
    }

    .pill-button {
      flex: 1;
      padding: 0.7rem 0.8rem;
      font-size: 0.85rem;
      justify-content: center;
      min-width: 0;
    }

    .generate,
    .upload {
      padding: 0.7rem;
      flex: 0 0 auto;
      width: 44px;
    }

    .pill-button span {
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }
  }
</style>
