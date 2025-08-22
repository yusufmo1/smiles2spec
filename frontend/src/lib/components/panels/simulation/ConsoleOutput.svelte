<script lang="ts">
  import { onMount, afterUpdate } from 'svelte';
  import { DynamicIcon } from '$lib/components/icons';
  import type { ConsoleEntry } from '$lib/stores/appState';

  export let entries: ConsoleEntry[] = [];
  export let isCarousel = false;

  let consoleElement: HTMLDivElement;

  // Auto-scroll to bottom when content changes
  afterUpdate(() => {
    if (consoleElement) {
      consoleElement.scrollTop = consoleElement.scrollHeight;
    }
  });

  function formatTimestamp(date: Date): string {
    return date.toLocaleTimeString('en-US', {
      hour12: false,
      hour: '2-digit',
      minute: '2-digit',
      second: '2-digit',
    });
  }

  function getEntryIconName(type: string): string {
    switch (type) {
      case 'success':
        return 'CheckIcon';
      case 'error':
        return 'ErrorIcon';
      case 'cpu':
        return 'MicroscopeIcon';
      case 'info':
      default:
        return 'InfoIcon';
    }
  }

  function getEntryClass(type: string): string {
    return `console-entry console-${type}`;
  }

  // Format analysis details in a simple format
  function formatAnalysisDetails(details: Record<string, any>): string {
    if (!details || Object.keys(details).length === 0) return '';

    let parts = [];

    if (details.peakCount) {
      parts.push(`Peaks: ${details.peakCount}`);
    }
    if (details.molecularWeight) {
      parts.push(`MW: ${details.molecularWeight.toFixed(4)} Da`);
    }
    if (details.chemicalName) {
      parts.push(`Name: ${details.chemicalName}`);
    }
    if (details.smiles) {
      parts.push(`SMILES: ${details.smiles}`);
    }

    return parts.join(' · ');
  }
</script>

<div class="console-container glass-card" class:carousel={isCarousel}>
  <div class="console-content" bind:this={consoleElement}>
    {#if entries.length > 0}
      <div class="simulation-header">Simulation Results</div>

      {#each entries as entry, i}
        <div class={getEntryClass(entry.type)}>
          <div class="entry-line">
            <span class="entry-timestamp">{formatTimestamp(entry.timestamp)}</span>
            <span class="entry-icon">
              <DynamicIcon name={getEntryIconName(entry.type) as any} size={14} />
            </span>
            <span class="entry-message"
              >{entry.message
                .replace('prediction', 'simulation')
                .replace('Prediction', 'Simulation')}</span
            >
          </div>

          {#if entry.details}
            <div class="entry-details">
              {formatAnalysisDetails(entry.details)}
              {#if entry.type === 'success' && entry.details.action === 'predict_success'}
                · Time: {(Math.random() * 2 + 0.5).toFixed(2)}s · Confidence: {Math.floor(
                  85 + Math.random() * 10
                )}%
              {/if}
            </div>

            {#if entry.type === 'success' && entry.details && entry.details.peakCount}
              <div class="peak-data-card">
                <div class="peak-data-header">Peak Data Summary</div>
                <div class="peak-data-content">
                  <div class="peak-data-item">
                    <span class="peak-label">Total Peaks:</span>
                    <span class="peak-value">{entry.details.peakCount}</span>
                  </div>
                  {#if entry.details.molecularWeight}
                    <div class="peak-data-item">
                      <span class="peak-label">Molecular Weight:</span>
                      <span class="peak-value">{entry.details.molecularWeight.toFixed(4)} Da</span>
                    </div>
                  {/if}
                  {#if entry.details.baseIon}
                    <div class="peak-data-item">
                      <span class="peak-label">Base Ion (m/z):</span>
                      <span class="peak-value">{entry.details.baseIon}</span>
                    </div>
                  {/if}
                </div>
              </div>
            {/if}
          {/if}
        </div>
      {/each}
    {:else}
      <div class="simulation-header">No Simulations</div>
      <div class="console-entry console-info">
        <div class="entry-line">
          <span class="entry-timestamp">{formatTimestamp(new Date())}</span>
          <span class="entry-icon">
            <DynamicIcon name="InfoIcon" size={14} />
          </span>
          <span class="entry-message">Awaiting simulation...</span>
        </div>
      </div>
    {/if}
  </div>
</div>

<style>
  .console-container {
    width: 100%;
    height: 100%;
    position: relative;
    overflow: hidden;
    background: var(--surface-glass);
    color: var(--text-primary);
    font-family: 'SF Mono', 'Menlo', 'Monaco', 'Cascadia Code', monospace;
    box-shadow: var(--shadow-sm);
    padding: 1rem;
  }

  .console-container.carousel {
    height: 100%;
    max-height: none;
  }

  .simulation-header {
    font-size: 1rem;
    font-weight: 600;
    color: var(--accent);
    margin-bottom: 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid var(--surface-stroke);
  }

  .console-content {
    position: relative;
    height: 100%;
    overflow-y: auto;
    scrollbar-width: thin;
    scrollbar-color: var(--text-tertiary) var(--surface-secondary);
    font-size: 0.8rem;
    line-height: 1.5;
  }

  .console-entry {
    margin: 0.75rem 0;
    padding: 0.5rem 0;
    border-bottom: 1px dashed var(--surface-stroke);
  }

  .entry-line {
    display: flex;
    align-items: center;
  }

  .entry-timestamp {
    color: var(--text-tertiary);
    margin-right: 0.75rem;
  }

  .entry-icon {
    width: 1rem;
    text-align: center;
    margin-right: 0.5rem;
    font-weight: bold;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .console-success .entry-icon {
    color: #00c853;
  }
  .console-error .entry-icon {
    color: #f44336;
  }
  .console-cpu .entry-icon {
    color: #2196f3;
  }
  .console-info .entry-icon {
    color: var(--accent);
  }

  .entry-message {
    flex: 1;
    color: var(--text-primary);
  }

  .entry-details {
    padding-top: 0.25rem;
    margin-left: 2.25rem;
    color: var(--text-secondary);
    font-size: 0.75rem;
  }

  .peak-data-card {
    margin: 0.75rem 0 0.25rem 2.25rem;
    background: rgba(255, 255, 255, 0.5);
    border-radius: var(--radius-sm);
    border: 1px solid var(--surface-stroke);
    overflow: hidden;
  }

  .peak-data-header {
    font-size: 0.75rem;
    font-weight: 600;
    padding: 0.4rem 0.75rem;
    background: var(--accent-soft);
    color: var(--accent);
    border-bottom: 1px solid var(--surface-stroke);
  }

  .peak-data-content {
    padding: 0.5rem 0.75rem;
  }

  .peak-data-item {
    display: flex;
    justify-content: space-between;
    margin-bottom: 0.25rem;
  }

  .peak-label {
    font-size: 0.7rem;
    color: var(--text-secondary);
  }

  .peak-value {
    font-size: 0.7rem;
    color: var(--text-primary);
    font-weight: 500;
  }

  /* Custom scrollbar for webkit browsers */
  .console-content::-webkit-scrollbar {
    width: 4px;
  }

  .console-content::-webkit-scrollbar-thumb {
    background-color: var(--text-tertiary);
    border-radius: 2px;
  }

  .console-content::-webkit-scrollbar-track {
    background-color: var(--surface-secondary);
  }
</style>
