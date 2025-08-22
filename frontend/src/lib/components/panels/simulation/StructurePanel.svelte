<script lang="ts">
  import { PlaceholderState } from '../components';
  import { DynamicIcon } from '$lib/components/icons';
  import { appState } from '$lib/stores/appState';
  import { smilesInputStore } from '$lib/components/smiles-input/stores/smilesInputStore';

  export let png: string | null = null;
  export let isCarousel = false;
  export let smiles: string = '';

  function loadToInput() {
    if (smiles) {
      smilesInputStore.setSmiles(smiles);
      // Add console feedback
      appState.addConsoleEntry({
        type: 'info',
        message: `SMILES loaded to input: ${smiles}`,
        details: { smiles, action: 'load_to_input' },
      });
    }
  }
</script>

<div class="structure-wrapper" class:carousel={isCarousel}>
  {#if png}
    <div class="structure-container">
      <img src={`data:image/png;base64,${png}`} alt="Molecular structure" />
      {#if smiles}
        <button class="load-button" on:click={loadToInput} title="Load SMILES to input">
          <DynamicIcon name="RefreshIcon" size={16} color="white" />
          <span>Load to Input</span>
        </button>
      {/if}
    </div>
  {:else}
    <PlaceholderState icon="MoleculeIcon" message="Awaiting molecule" />
  {/if}
</div>

<style>
  .structure-wrapper {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 100%;
    height: 100%;
    overflow: hidden;
  }

  .structure-wrapper.carousel {
    height: 100%;
    max-height: none;
  }

  .structure-container {
    position: relative;
    display: flex;
    align-items: center;
    justify-content: center;
    width: 100%;
    height: 100%;
  }

  img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
  }

  .load-button {
    position: absolute;
    top: 12px;
    right: 12px;
    background: var(--accent);
    color: white;
    border: none;
    border-radius: 20px;
    padding: 8px 12px;
    font-size: 0.8rem;
    font-weight: 500;
    display: flex;
    align-items: center;
    gap: 6px;
    cursor: pointer;
    transition: all 0.2s ease;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
    opacity: 0.9;
  }

  .load-button:hover {
    background: var(--accent-secondary);
    transform: translateY(-1px);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
    opacity: 1;
  }

  .load-button:active {
    transform: translateY(0);
  }
</style>
