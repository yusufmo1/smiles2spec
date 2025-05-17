<script>
  import { createEventDispatcher } from 'svelte';
  import { exportMsp, exportMspBatch } from '../services/api.js';
  
  export let name = "";
  export let index = 0;   // 1-based
  export let total = 0;
  export let smilesList = [];
  export let currentSmiles = "";
  
  const dispatch = createEventDispatcher();
  
  async function exportMspData() {
    try {
      let blob;
      
      if (total > 1) {
        // Export all entries as batch
        blob = await exportMspBatch(smilesList);
      } else {
        // Export single entry
        blob = await exportMsp(currentSmiles);
      }
      
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = total > 1 ? `batch_export.msp` : `${currentSmiles}.msp`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      alert(`Could not export MSP: ${e.message}`);
    }
  }
</script>

<div class="navigator glass-card">
  <div class="spacer"></div>
  
  <div class="center-group">
    <button class="arrow" on:click={() => dispatch('prev')} disabled={total<=1}>◀</button>
    <div class="name-container">
      <span class="name">{name || "—"}</span>
      {#if total}
        <span class="count">{index}/{total}</span>
      {/if}
    </div>
    <button class="arrow" on:click={() => dispatch('next')} disabled={total<=1}>▶</button>
  </div>
  
  <div class="button-container">
    <button class="export-btn" 
            on:click={exportMspData} 
            disabled={!currentSmiles} 
            title={total > 1 ? "Export all entries as MSP" : "Export current entry as MSP"}>
      <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
        <path d="M4 18h16v2H4v-2zM4 14h16v2H4v-2zM8 4v6h3v6h2v-6h3V4H8z"/>
      </svg>
      {total > 1 ? 'Export All' : 'Export MSP'}
    </button>
  </div>
</div>

<style>
.navigator {
  margin-top: -0.5rem;
  width: 100%;
  padding: 0.75rem 1.25rem;
  display: grid;
  grid-template-columns: 1fr auto 1fr;
  align-items: center;
  font-weight: 500;
  font-size: 0.95rem;
  border-radius: var(--enforce-pill);
}

.center-group {
  display: flex;
  align-items: center;
  gap: 0.75rem;
}

.spacer {
  /* Empty spacer for grid layout */
}

.button-container {
  display: flex;
  justify-content: flex-end;
}

.name-container {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  text-align: center;
}

.name {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  max-width: 250px;
}

.count { 
  color: var(--text-tertiary); 
  font-size: 0.8rem; 
}

.arrow {
  background: transparent;
  border: none;
  font-size: 1.1rem;
  color: var(--accent);
  cursor: pointer;
  padding: 0;
  width: 1.5rem;
  height: 1.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
}

.arrow:disabled { 
  opacity: 0.25; 
  cursor: default; 
}

.export-btn {
  background: var(--accent-soft);
  border: none;
  color: var(--accent);
  font-size: 0.8rem;
  font-weight: 500;
  padding: 0.35rem 0.75rem;
  border-radius: var(--enforce-pill);
  cursor: pointer;
  display: flex;
  align-items: center;
  gap: 0.35rem;
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