<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  let activeTab = 'architecture';
  let progress = 0;

  import { onMount } from 'svelte';

  onMount(() => {
    const interval = setInterval(() => {
      progress = (progress + 1) % 101;
    }, 30);

    return () => clearInterval(interval);
  });
</script>

<div class="ml-engine-panel">
  <div class="tabs">
    <button
      class="tab"
      class:active={activeTab === 'architecture'}
      on:click={() => (activeTab = 'architecture')}
    >
      Architecture
    </button>
    <button
      class="tab"
      class:active={activeTab === 'training'}
      on:click={() => (activeTab = 'training')}
    >
      Training
    </button>
  </div>

  {#if activeTab === 'architecture'}
    <div class="architecture-view">
      <div class="model-visualization">
        <div class="layer input-layer">
          <div class="layer-icon"><DynamicIcon name="MoleculeIcon" size={24} /></div>
          <span>Input Layer</span>
          <div class="layer-info">8000+ features</div>
        </div>

        <div class="connection">
          <div class="flow-line"></div>
        </div>

        <div class="layer hidden-layer">
          <div class="layer-icon"><DynamicIcon name="CpuIcon" size={24} /></div>
          <span>DRF Engine</span>
          <div class="layer-info">Differentiable Random Forest</div>
        </div>

        <div class="connection">
          <div class="flow-line"></div>
        </div>

        <div class="layer output-layer">
          <div class="layer-icon"><DynamicIcon name="SpectrumIcon" size={24} /></div>
          <span>Output Layer</span>
          <div class="layer-info">Mass spectrum</div>
        </div>
      </div>

      <div class="model-stats">
        <div class="stat-card">
          <DynamicIcon name="FlashIcon" size={20} color="var(--accent)" />
          <div class="stat-value">&gt;100ms</div>
          <div class="stat-label">Inference Time</div>
        </div>
        <div class="stat-card">
          <DynamicIcon name="TargetIcon" size={20} color="var(--accent)" />
          <div class="stat-value">SOTA</div>
          <div class="stat-label">Performance</div>
        </div>
      </div>
    </div>
  {:else}
    <div class="training-view">
      <div class="training-progress">
        <div class="progress-circle">
          <svg viewBox="0 0 100 100">
            <circle
              cx="50"
              cy="50"
              r="45"
              fill="none"
              stroke="rgba(120, 121, 255, 0.1)"
              stroke-width="5"
            />
            <circle
              cx="50"
              cy="50"
              r="45"
              fill="none"
              stroke="var(--accent)"
              stroke-width="5"
              stroke-dasharray={`${progress * 2.83} 283`}
              stroke-linecap="round"
              transform="rotate(-90 50 50)"
            />
          </svg>
          <div class="progress-value">{progress}%</div>
        </div>
      </div>

      <div class="training-info">
        <div class="info-item">
          <DynamicIcon name="DatabaseIcon" size={20} />
          <span>50,000+ Spectra</span>
        </div>
        <div class="info-item">
          <DynamicIcon name="CpuIcon" size={20} />
          <span>Teacher-Student Framework</span>
        </div>
      </div>
    </div>
  {/if}
</div>

<style>
  .ml-engine-panel {
    height: 100%;
    display: flex;
    flex-direction: column;
  }

  .tabs {
    display: flex;
    gap: 0.5rem;
    padding: 1rem;
    justify-content: center;
  }

  .tab {
    padding: 0.5rem 1.5rem;
    background: rgba(255, 255, 255, 0.6);
    border: 1px solid rgba(120, 121, 255, 0.3);
    border-radius: 20px;
    color: var(--text-primary);
    font-size: 0.9rem;
    cursor: pointer;
    transition: all 0.2s ease;
  }

  .tab:hover {
    background: rgba(120, 121, 255, 0.1);
  }

  .tab.active {
    background: var(--accent);
    color: white;
    border-color: var(--accent);
  }

  .architecture-view {
    flex: 1;
    display: flex;
    flex-direction: column;
    gap: 2rem;
    padding: 2rem;
  }

  .model-visualization {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 1rem;
  }

  .layer {
    background: rgba(255, 255, 255, 0.8);
    border-radius: 16px;
    padding: 1.5rem 1rem;
    text-align: center;
    min-width: 100px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0.5rem;
  }

  .layer-icon {
    color: var(--accent);
  }

  .layer span {
    font-weight: 600;
    font-size: 0.9rem;
  }

  .layer-info {
    font-size: 0.75rem;
    color: var(--text-secondary);
  }

  .connection {
    width: 40px;
    height: 2px;
    position: relative;
  }

  .flow-line {
    width: 100%;
    height: 100%;
    background: linear-gradient(90deg, var(--accent), var(--accent-secondary));
    position: relative;
    overflow: hidden;
  }

  .flow-line::after {
    content: '';
    position: absolute;
    top: -8px;
    left: -20px;
    width: 20px;
    height: 20px;
    background: white;
    border-radius: 50%;
    box-shadow: 0 0 10px var(--accent);
    animation: flow 2s linear infinite;
  }

  @keyframes flow {
    to {
      left: 100%;
    }
  }

  .model-stats {
    display: flex;
    gap: 1rem;
    justify-content: center;
  }

  .stat-card {
    background: rgba(120, 121, 255, 0.1);
    border-radius: 12px;
    padding: 1rem;
    text-align: center;
    min-width: 100px;
  }

  .stat-value {
    font-size: 1.5rem;
    font-weight: 700;
    color: var(--accent);
    margin: 0.5rem 0;
  }

  .stat-label {
    font-size: 0.8rem;
    color: var(--text-secondary);
  }

  .training-view {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    gap: 2rem;
    padding: 2rem;
  }

  .training-progress {
    position: relative;
    width: 150px;
    height: 150px;
  }

  .progress-circle svg {
    width: 100%;
    height: 100%;
  }

  .progress-value {
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    font-size: 2rem;
    font-weight: 700;
    color: var(--accent);
  }

  .training-info {
    display: flex;
    gap: 2rem;
  }

  .info-item {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    color: var(--text-secondary);
    font-size: 0.9rem;
  }
</style>
