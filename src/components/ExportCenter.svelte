<script>
  import { exportMsp, exportMspBatch } from '../services/api.js';
  
  // Export center component for all export functionality
  export let spectrumData = null;
  export let peakData = [];
  export let structurePNG = null;
  export let smiles = "";
  export let chemicalName = "";
  export let smilesList = [];
  
  // Track available export options
  $: hasSpectrum = spectrumData !== null;
  $: hasPeaks = Array.isArray(peakData) && peakData.length > 0;
  $: hasStructure = structurePNG !== null;
  $: hasChemicalInfo = !!smiles || !!chemicalName;
  $: hasBatchData = Array.isArray(smilesList) && smilesList.length > 1;
  
  async function exportMspData() {
    try {
      let blob;
      
      if (hasBatchData) {
        // Export all entries as batch
        blob = await exportMspBatch(smilesList);
      } else if (smiles) {
        // Export single entry
        blob = await exportMsp(smiles);
      } else {
        throw new Error("No chemical data available to export");
      }
      
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = hasBatchData ? `batch_export.msp` : `${smiles || chemicalName || "export"}.msp`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      alert(`Could not export MSP: ${e.message}`);
    }
  }
</script>

<div class="export-center">
  <div class="export-status">
    {#if hasChemicalInfo}
      <span class="status-item available">✓ Chemical data</span>
    {:else}
      <span class="status-item unavailable">✗ No chemical data</span>
    {/if}
    
    {#if hasSpectrum}
      <span class="status-item available">✓ Spectrum data</span>
    {:else}
      <span class="status-item unavailable">✗ No spectrum data</span>
    {/if}
    
    {#if hasStructure}
      <span class="status-item available">✓ Structure image</span>
    {:else}
      <span class="status-item unavailable">✗ No structure image</span>
    {/if}
  </div>
  
  <div class="export-options">
    <button class="export-btn" 
            on:click={exportMspData} 
            disabled={!smiles && (!smilesList || smilesList.length === 0)} 
            title="Export as MSP format">
      <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
        <path d="M4 18h16v2H4v-2zM4 14h16v2H4v-2zM8 4v6h3v6h2v-6h3V4H8z"/>
      </svg>
      {hasBatchData ? 'Export All MSP' : 'Export MSP'}
    </button>
  </div>
</div>

<style>
  .export-center {
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
    gap: 1.5rem;
  }
  
  .export-status {
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    justify-content: center;
    font-size: 0.75rem;
  }
  
  .status-item {
    padding: 0.35rem 0.75rem;
    border-radius: var(--enforce-pill);
  }
  
  .available {
    background: rgba(52, 199, 89, 0.15);
    color: #34c759;
  }
  
  .unavailable {
    background: rgba(142, 142, 147, 0.15);
    color: #8e8e93;
  }
  
  .export-options {
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
    justify-content: center;
  }
  
  .export-btn {
    background: var(--accent-soft);
    border: none;
    color: var(--accent);
    font-size: 0.8rem;
    font-weight: 500;
    padding: 0.5rem 1rem;
    border-radius: var(--enforce-pill);
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    transition: all 0.2s;
  }

  .export-btn:hover:not(:disabled) {
    background: var(--accent);
    color: white;
  }

  .export-btn:disabled {
    opacity: 0.4;
    cursor: not-allowed;
  }
</style> 