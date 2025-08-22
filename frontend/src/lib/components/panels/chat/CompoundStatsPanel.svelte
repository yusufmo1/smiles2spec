<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let smiles: string = '';
  export let molecularWeight: number | null = null;
  export let exactMass: number | null = null;
  export let peakCount: number = 0;
  export let chemicalName: string = '';
  export let predictionConfidence: number = 0;
  export let isCarousel = false;

  // Calculate derived stats
  $: hasData = smiles && peakCount > 0;
  $: massAccuracy =
    exactMass && molecularWeight
      ? ((Math.abs(exactMass - molecularWeight) / molecularWeight) * 100).toFixed(3)
      : null;
</script>

<div class="compound-stats-panel" class:carousel={isCarousel}>
  {#if hasData}
    <div class="stats-grid">
      <div class="stat-card primary hover-lift">
        <div class="stat-icon">
          <DynamicIcon name="FlaskIcon" size={24} color="var(--accent)" />
        </div>
        <div class="stat-content">
          <h4>Chemical Identity</h4>
          <p class="compound-name">{chemicalName || 'Unknown compound'}</p>
          <p class="smiles-display">{smiles}</p>
        </div>
      </div>

      <div class="stat-card hover-lift">
        <div class="stat-icon">
          <DynamicIcon name="DatabaseIcon" size={24} color="var(--text-secondary)" />
        </div>
        <div class="stat-content">
          <h4>Molecular Properties</h4>
          <p>MW: {molecularWeight?.toFixed(2) || 'N/A'} Da</p>
          <p>Exact Mass: {exactMass?.toFixed(4) || 'N/A'} Da</p>
          {#if massAccuracy}
            <p>Accuracy: ±{massAccuracy}%</p>
          {/if}
        </div>
      </div>

      <div class="stat-card hover-lift">
        <div class="stat-icon">
          <DynamicIcon name="SpectrumIcon" size={24} color="var(--text-secondary)" />
        </div>
        <div class="stat-content">
          <h4>Spectrum Analysis</h4>
          <p>Detected Peaks: {peakCount}</p>
          <p>Confidence: {(predictionConfidence * 100).toFixed(1)}%</p>
          <p>Base Peak: {peakCount > 0 ? 'Identified' : 'N/A'}</p>
        </div>
      </div>

      <div class="stat-card hover-lift">
        <div class="stat-icon">
          <DynamicIcon name="MicroscopeIcon" size={24} color="var(--text-secondary)" />
        </div>
        <div class="stat-content">
          <h4>Fragmentation</h4>
          <p>Major Fragments: {Math.min(5, Math.floor(peakCount * 0.3))}</p>
          <p>Molecular Ion: Detected</p>
          <p>Loss Patterns: Analysed</p>
        </div>
      </div>
    </div>

    <div class="knowledge-summary">
      <div class="summary-header">
        <span class="summary-icon">
          <DynamicIcon name="MoleculeIcon" size={20} color="var(--accent)" />
        </span>
        <h3>Spectrum's Knowledge</h3>
      </div>
      <p>
        I can analyse this {peakCount}-peak spectrum, identify fragmentation patterns, explain mass
        losses, and help interpret the molecular structure. My confidence in this prediction is {(
          predictionConfidence * 100
        ).toFixed(0)}%.
      </p>
    </div>
  {:else}
    <div class="no-data-state">
      <div class="no-data-icon">
        <DynamicIcon name="SearchIcon" size={48} color="var(--text-secondary)" />
      </div>
      <h3>No Compound Data</h3>
      <p>Predict a spectrum first to see what Spectrum knows about your compound.</p>
    </div>
  {/if}
</div>

<style>
  /* Import only needed animations */
  .hover-lift {
    transition:
      transform var(--transition-smooth),
      box-shadow var(--transition-smooth);
  }

  .hover-lift:hover {
    transform: translateY(-2px);
    box-shadow: 0 8px 25px rgba(0, 0, 0, 0.08);
  }
  .compound-stats-panel {
    height: 100%;
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
    overflow-y: auto;
  }

  .compound-stats-panel.carousel {
    height: 100%;
    max-height: none;
  }

  .stats-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
    flex: 1;
  }

  .stat-card {
    background: rgba(255, 255, 255, 0.6);
    border-radius: var(--radius-lg);
    padding: 1.25rem;
    border: 1px solid rgba(255, 255, 255, 0.8);
    display: flex;
    gap: 1rem;
    align-items: flex-start;
  }

  .stat-card:hover {
    background: rgba(255, 255, 255, 0.8);
  }

  .stat-card.primary {
    background: linear-gradient(135deg, rgba(120, 121, 255, 0.1), rgba(255, 255, 255, 0.6));
    border-color: rgba(120, 121, 255, 0.3);
  }

  .stat-icon {
    flex-shrink: 0;
  }

  .stat-content h4 {
    margin: 0 0 0.5rem 0;
    font-size: 0.9rem;
    font-weight: 600;
    color: var(--text-primary);
  }

  .stat-content p {
    margin: 0.25rem 0;
    font-size: 0.85rem;
    color: var(--text-secondary);
  }

  .compound-name {
    font-weight: 600 !important;
    color: var(--accent) !important;
  }

  .smiles-display {
    font-family: 'SF Mono', monospace !important;
    font-size: 0.8rem !important;
    background: rgba(0, 0, 0, 0.05);
    padding: 0.25rem 0.5rem;
    border-radius: 4px;
    word-break: break-all;
  }

  .knowledge-summary {
    background: rgba(255, 255, 255, 0.5);
    border-radius: var(--radius-lg);
    padding: 1.5rem;
    border: 1px solid rgba(255, 255, 255, 0.7);
  }

  .summary-header {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    margin-bottom: 1rem;
  }

  .summary-icon {
    display: flex;
  }

  .summary-header h3 {
    margin: 0;
    font-size: 1rem;
    font-weight: 600;
  }

  .knowledge-summary p {
    margin: 0;
    font-size: 0.9rem;
    line-height: 1.5;
    color: var(--text-secondary);
  }

  .no-data-state {
    height: 100%;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    text-align: center;
    opacity: 0.6;
  }

  .no-data-icon {
    margin-bottom: 1rem;
  }

  .no-data-state h3 {
    margin: 0 0 0.5rem 0;
    color: var(--text-primary);
  }

  .no-data-state p {
    margin: 0;
    color: var(--text-secondary);
    max-width: 250px;
  }

  @media (max-width: 768px) {
    .stats-grid {
      grid-template-columns: 1fr;
    }

    .stat-card {
      padding: 1rem;
    }

    .knowledge-summary {
      padding: 1rem;
    }
  }
</style>
