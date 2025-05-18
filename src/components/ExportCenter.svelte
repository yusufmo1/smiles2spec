<script>
  // Export center component - placeholder for future functionality
  export let spectrumData = null;
  export let peakData = [];
  export let structurePNG = null;
  export let smiles = "";
  export let chemicalName = "";
</script>

<div class="export-center">
  <div class="placeholder-message">
    <p>Export functionality will be implemented in a future update.</p>
  </div>
</div>

<style>
  .export-center {
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
  }
  
  .placeholder-message {
    color: var(--text-secondary);
    font-style: italic;
    text-align: center;
  }
</style> 