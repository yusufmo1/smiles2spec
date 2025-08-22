<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let spectrumData: any = null;
  export let peakData: Array<{ mz: number; intensity: number }> = [];
  export let structurePNG: string | null = null;
  export let smiles: string = '';
  export let chemicalName: string = '';
  export let smilesList: string[] = [];
  export let isCarousel = false;

  // Track available export options
  $: hasSpectrum = spectrumData !== null;
  $: hasPeaks = Array.isArray(peakData) && peakData.length > 0;
  $: hasStructure = structurePNG !== null;
  $: hasChemicalInfo = !!smiles || !!chemicalName;
  $: hasBatchData = Array.isArray(smilesList) && smilesList.length > 1;

  function exportJSON() {
    try {
      const exportData = {
        chemical_info: {
          smiles,
          chemical_name: chemicalName,
        },
        spectrum_data: spectrumData,
        peak_data: peakData,
        has_structure: hasStructure,
        export_timestamp: new Date().toISOString(),
      };

      const blob = new Blob([JSON.stringify(exportData, null, 2)], {
        type: 'application/json',
      });

      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${smiles || chemicalName || 'export'}_spectrum_data.json`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      alert(`Could not export JSON: ${(e as Error).message}`);
    }
  }

  function exportCSV() {
    try {
      if (!hasPeaks) {
        throw new Error('No peak data available to export');
      }

      let csvContent = 'M/Z,Intensity,Relative_Intensity\n';
      const maxIntensity = Math.max(...peakData.map((p) => p.intensity));

      peakData.forEach((peak) => {
        const relativeIntensity = ((peak.intensity / maxIntensity) * 100).toFixed(2);
        csvContent += `${peak.mz.toFixed(4)},${peak.intensity.toFixed(6)},${relativeIntensity}\n`;
      });

      const blob = new Blob([csvContent], { type: 'text/csv' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${smiles || chemicalName || 'export'}_peak_data.csv`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      alert(`Could not export CSV: ${(e as Error).message}`);
    }
  }

  function exportMSP() {
    try {
      if (!hasChemicalInfo || !hasPeaks) {
        throw new Error('Need both chemical info and peak data for MSP export');
      }

      let mspContent = `NAME: ${chemicalName || smiles}\n`;
      mspContent += `SMILES: ${smiles}\n`;
      mspContent += `MW: Unknown\n`;
      mspContent += `Num Peaks: ${peakData.length}\n\n`;

      peakData.forEach((peak) => {
        mspContent += `${peak.mz.toFixed(4)}\t${peak.intensity.toFixed(6)}\n`;
      });

      const blob = new Blob([mspContent], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${smiles || chemicalName || 'export'}.msp`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      alert(`Could not export MSP: ${(e as Error).message}`);
    }
  }
</script>

<div class="export-center" class:carousel={isCarousel}>
  <div class="export-status">
    {#if hasChemicalInfo}
      <span class="status-item available">
        <span class="status-icon">
          <DynamicIcon name="CheckIcon" size={12} />
        </span> Chemical data
      </span>
    {:else}
      <span class="status-item unavailable">
        <span class="status-icon">
          <DynamicIcon name="InfoIcon" size={12} />
        </span> No chemical data
      </span>
    {/if}

    {#if hasSpectrum}
      <span class="status-item available">
        <span class="status-icon">
          <DynamicIcon name="CheckIcon" size={12} />
        </span> Spectrum data
      </span>
    {:else}
      <span class="status-item unavailable">
        <span class="status-icon">
          <DynamicIcon name="InfoIcon" size={12} />
        </span> No spectrum data
      </span>
    {/if}

    {#if hasPeaks}
      <span class="status-item available">
        <span class="status-icon">
          <DynamicIcon name="CheckIcon" size={12} />
        </span> Peak data
      </span>
    {:else}
      <span class="status-item unavailable">
        <span class="status-icon">
          <DynamicIcon name="InfoIcon" size={12} />
        </span> No peak data
      </span>
    {/if}
  </div>

  <div class="export-options">
    <button
      class="export-btn hover-lift-subtle"
      on:click={exportJSON}
      disabled={!hasChemicalInfo && !hasSpectrum}
      title="Export raw data as JSON"
    >
      <DynamicIcon name="ExportIcon" size={16} />
      Export JSON
    </button>

    <button
      class="export-btn hover-lift-subtle"
      on:click={exportCSV}
      disabled={!hasPeaks}
      title="Export peak data as CSV"
    >
      <DynamicIcon name="ExportIcon" size={16} />
      Export CSV
    </button>

    <button
      class="export-btn hover-lift-subtle"
      on:click={exportMSP}
      disabled={!hasChemicalInfo || !hasPeaks}
      title="Export as MSP format"
    >
      <DynamicIcon name="ExportIcon" size={16} />
      Export MSP
    </button>
  </div>
</div>

<style>
  /* Import only needed animations */
  .hover-lift-subtle {
    transition:
      transform var(--transition-smooth),
      box-shadow var(--transition-smooth);
  }

  .hover-lift-subtle:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.06);
  }
  .export-center {
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
    gap: 1.5rem;
    padding: 1rem;
  }

  .export-center.carousel {
    height: 100%;
    max-height: none;
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
    border-radius: var(--radius-pill);
    display: flex;
    align-items: center;
    gap: 0.25rem;
  }

  .status-icon {
    flex-shrink: 0;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .available {
    background: rgba(52, 199, 89, 0.15);
    color: #34c759;
    border: 1px solid rgba(52, 199, 89, 0.25);
  }

  .unavailable {
    background: rgba(142, 142, 147, 0.15);
    color: #8e8e93;
    border: 1px solid rgba(142, 142, 147, 0.25);
  }

  .export-options {
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
    justify-content: center;
    width: 100%;
  }

  .export-btn {
    background: var(--accent-soft);
    border: 1px solid rgba(120, 121, 255, 0.25);
    color: var(--accent);
    font-size: 0.8rem;
    font-weight: 500;
    padding: 0.65rem 1.25rem;
    border-radius: var(--radius-pill);
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    box-shadow: 0 2px 8px rgba(120, 121, 255, 0.1);
  }

  .export-btn:hover:not(:disabled) {
    background: var(--accent);
    color: white;
    border-color: var(--accent);
  }

  .export-btn:disabled {
    opacity: 0.4;
    cursor: not-allowed;
    box-shadow: none;
  }

  .export-btn :global(svg) {
    flex-shrink: 0;
  }
</style>
