<script lang="ts">
  import { panelStore, simulationPanels, appState } from '$lib/stores/index';
  import Panel from '../Panel.svelte';
  import type { Panel as PanelType } from '$lib/stores/index';
  import { devLog } from '$lib/utils/performance';

  // Accept panels prop or use default simulation panels
  export let panels: PanelType[] | null = null;
  export let currentPage: string = 'home';

  // Data props for simulation panels
  export let spectrumData: any = null;
  export let peakData: any[] = [];
  export let structurePNG: string | null = null;
  export let currentSmiles: string = '';
  export let currentName: string = '';
  export let consoleText: string = '';
  export let hasFirstPrediction: boolean = false;
  export let smilesList: string[] = [];

  // Log when props change
  $: {
    devLog.log(`PanelGrid props updated:`, {
      spectrumData: !!spectrumData,
      peakData: peakData.length,
      structurePNG: !!structurePNG,
      currentSmiles,
      currentName,
      hasFirstPrediction,
      smilesList: smilesList.length,
    });
  }

  // Update panel props reactively when data changes
  $: if (currentPage === 'home' || currentPage === 'spectral-simulation') {
    devLog.log('Updating panel props for simulation page');

    panelStore.updatePanelProps('spectrum', {
      spectrumData,
      isCarousel: false,
    });

    panelStore.updatePanelProps('fragmentation', {
      data: spectrumData,
      isCarousel: false,
    });

    panelStore.updatePanelProps('structure', {
      png: structurePNG,
      smiles: currentSmiles,
      isCarousel: false,
    });

    panelStore.updatePanelProps('peaks', {
      peaks: peakData,
      smiles: currentSmiles,
      isCarousel: false,
    });

    panelStore.updatePanelProps('console', {
      output: consoleText,
      entries: $appState.consoleEntries, // Pass detailed entries
      isCarousel: false,
    });

    panelStore.updatePanelProps('export', {
      spectrumData,
      peakData,
      structurePNG,
      smiles: currentSmiles,
      chemicalName: currentName,
      smilesList,
      isCarousel: false,
    });

    panelStore.updatePanelProps('chat', {
      hasSmilesPrediction: hasFirstPrediction,
      currentSmiles,
      isCarousel: false,
    });
  }

  // Use provided panels or default simulation panels
  $: currentPanels = panels || $simulationPanels;

  // Log current panels
  $: {
    devLog.log(
      'PanelGrid currentPanels:',
      currentPanels.length,
      currentPanels.map((p) => p.id)
    );
  }
</script>

<div class="panel-grid">
  {#if currentPanels && currentPanels.length > 0}
    {#each Array(Math.ceil(currentPanels.length / 2)) as _, rowIndex}
      <div class="row">
        {#if currentPanels[rowIndex * 2]}
          <div class="col-half">
            <Panel {...currentPanels[rowIndex * 2]} />
          </div>
        {/if}
        {#if currentPanels[rowIndex * 2 + 1]}
          <div class="col-half">
            <Panel {...currentPanels[rowIndex * 2 + 1]} />
          </div>
        {/if}
      </div>
    {/each}
  {:else}
    <!-- Fallback when no panels are available -->
    <div class="no-panels">
      <p>Loading panels...</p>
    </div>
  {/if}
</div>

<style>
  .panel-grid {
    display: flex;
    flex-direction: column;
    gap: 2rem;
    width: 100%;
    min-height: 200px; /* Ensure some height even when loading */
  }

  .row {
    display: flex;
    gap: 3rem;
    width: 100%;
  }

  .col-half {
    flex: 1;
    min-width: 0;
  }

  .no-panels {
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 200px;
    color: var(--text-secondary);
    font-style: italic;
  }

  @media (max-width: 1024px) {
    .row {
      flex-direction: column;
      gap: 2rem;
    }

    .panel-grid {
      gap: 1.5rem;
    }
  }

  @media (max-width: 768px) {
    .panel-grid {
      gap: 1rem;
    }

    .row {
      gap: 1rem;
    }
  }
</style>
