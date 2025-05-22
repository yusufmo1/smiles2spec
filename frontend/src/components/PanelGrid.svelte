<script>
  import { panelDefinitions, updatePanelProps } from '../stores/panelData.js';
  import Panel from './Panel.svelte';

  export let spectrumData = null;
  export let peakData = [];
  export let structurePNG = null;
  export let currentSmiles = "";
  export let currentName = "";
  export let consoleText = "";
  export let hasFirstPrediction = false;
  export let smilesList = [];

  $: {
    updatePanelProps('spectrum', { spectrumData });
    updatePanelProps('fragments', { data: spectrumData });
    updatePanelProps('structure', { png: structurePNG });
    updatePanelProps('peaks', { peaks: peakData, smiles: currentSmiles });
    updatePanelProps('console', { output: consoleText });
    updatePanelProps('chat', {
      hasSmilesPrediction: hasFirstPrediction,
      currentSmiles,
      smilesList
    });
  }
</script>

<div class="panel-grid">
  <div class="row">
    <div class="col-half">
      <Panel {...$panelDefinitions[0]} clickable={true} />
    </div>
    <div class="col-half">
      <Panel {...$panelDefinitions[1]} clickable={true} />
    </div>
  </div>

  <div class="row">
    <div class="col-half">
      <Panel {...$panelDefinitions[2]} clickable={true} />
    </div>
    <div class="col-half">
      <Panel {...$panelDefinitions[3]} clickable={true} />
    </div>
  </div>

  <div class="row">
    <div class="col-half">
      <Panel {...$panelDefinitions[4]} clickable={true} />
    </div>
    <div class="col-half">
      <Panel {...$panelDefinitions[5]} clickable={true} />
    </div>
  </div>
</div>

<style>
  .panel-grid {
    display: flex;
    flex-direction: column;
    gap: 2rem;
  }

  .row {
    display: flex;
    flex-wrap: wrap;
    margin: 0;
  }

  .col-half {
    width: 50%;
    padding: 0 1.5rem;
  }

  @media (max-width: 1024px) {
    .col-half {
      width: 100%;
      margin-bottom: 2rem;
    }
  }
</style>
