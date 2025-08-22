<script lang="ts">
  import { onMount, afterUpdate } from 'svelte';
  import { browser } from '$app/environment';
  import { DynamicIcon } from '$lib/components/icons';

  interface Peak {
    mz: number;
    intensity: number;
  }

  export let peaks: Peak[] = [];
  export let smiles: string = '';
  export let isCarousel = false;

  let sortField: 'mz' | 'intensity' = 'intensity';
  let sortDirection: 'asc' | 'desc' = 'desc';

  // Fragment description lookup based on common fragment patterns in mass spectrometry
  const fragmentDescriptions: Record<number, string> = {
    15: 'CH₃• (Methyl)',
    27: 'C₂H₃• (Vinyl)',
    28: 'CO/C₂H₄ (Carbonyl/Ethylene)',
    29: 'C₂H₅• (Ethyl)',
    30: 'CH₂NH₂• (Aminomethyl)',
    31: 'CH₂OH• (Hydroxymethyl)',
    39: 'C₃H₃• (Propynyl)',
    41: 'C₃H₅• (Allyl)',
    42: 'C₃H₆/CH₂CO (Propylene/Ketene)',
    43: 'C₃H₇•/CH₃CO• (Propyl/Acetyl)',
    44: 'CO₂/C₂H₄NH₂ (Carbon dioxide/Aminoethyl)',
    45: 'CH₃CHOH• (Ethanol fragment)',
    55: 'C₄H₇• (Butene fragment)',
    56: 'C₄H₈ (Butene)',
    57: 'C₄H₉• (Butyl)',
    69: 'C₅H₉• (Cyclopentenyl)',
    71: 'C₅H₁₁• (Pentyl)',
    77: 'C₆H₅• (Phenyl)',
    91: 'C₇H₇• (Tropylium)',
    105: 'C₈H₉• (Methyltropylium)',
  };

  // Custom fragment descriptions for specific SMILES patterns (could be expanded)
  function getCustomFragments(): Record<number, string> {
    if (!smiles) return {};

    const customFragments: Record<number, string> = {};

    // Check for common functional groups based on SMILES
    if (smiles.includes('O')) {
      // Compounds with oxygen might have these fragments
      customFragments[31] = 'CH₂OH• (Hydroxymethyl)';
      customFragments[45] = 'CH₃CHOH• (Ethanol fragment)';
    }

    if (smiles.includes('N')) {
      // Compounds with nitrogen might have these fragments
      customFragments[30] = 'CH₂NH₂• (Aminomethyl)';
      customFragments[44] = 'C₂H₄NH₂ (Aminoethyl)';
    }

    if (smiles.includes('c1ccccc1')) {
      // Compounds with benzene rings
      customFragments[77] = 'C₆H₅• (Phenyl)';
      customFragments[91] = 'C₇H₇• (Benzyl/Tropylium)';
    }

    return customFragments;
  }

  function getFragmentDescription(mz: number): string {
    const rounded = Math.round(mz);
    const customFrags = getCustomFragments();
    return customFrags[rounded] || fragmentDescriptions[rounded] || '';
  }

  function toggleSort(field: 'mz' | 'intensity') {
    if (sortField === field) {
      sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
    } else {
      sortField = field;
      sortDirection = 'desc';
    }
  }

  $: sortedPeaks = [...peaks].sort((a, b) => {
    const factor = sortDirection === 'asc' ? 1 : -1;
    return factor * (a[sortField] - b[sortField]);
  });

  $: baseIntensity = peaks.length > 0 ? Math.max(...peaks.map((p) => p.intensity)) : 1;
</script>

<div class="peak-table-container" class:carousel={isCarousel}>
  {#if peaks.length > 0}
    <div class="table-wrapper">
      <table class="peak-table">
        <thead>
          <tr>
            <th scope="col" on:click={() => toggleSort('mz')} class:active={sortField === 'mz'}>
              m/z
              <span class="sort-indicator"
                >{sortField === 'mz' ? (sortDirection === 'asc' ? '↑' : '↓') : ''}</span
              >
            </th>
            <th
              scope="col"
              on:click={() => toggleSort('intensity')}
              class:active={sortField === 'intensity'}
            >
              Intensity
              <span class="sort-indicator"
                >{sortField === 'intensity' ? (sortDirection === 'asc' ? '↑' : '↓') : ''}</span
              >
            </th>
            <th scope="col">Relative %</th>
            <th scope="col">Fragment</th>
          </tr>
        </thead>
        <tbody>
          {#each sortedPeaks as peak, i}
            <tr>
              <td>{peak.mz.toFixed(2)}</td>
              <td class="intensity-cell">
                <div
                  class="intensity-bar"
                  style="--percentage: {(peak.intensity / baseIntensity) * 100}%"
                >
                  {peak.intensity.toFixed(4)}
                </div>
              </td>
              <td>{((peak.intensity / baseIntensity) * 100).toFixed(1)}%</td>
              <td>{getFragmentDescription(peak.mz)}</td>
            </tr>
          {/each}
        </tbody>
      </table>
    </div>
  {:else}
    <div class="empty-state">
      <div class="empty-content">
        <DynamicIcon
          name="DataTableIcon"
          size={36}
          color="var(--text-secondary)"
          className="empty-icon"
        />
        <p>No peak data available</p>
      </div>
    </div>
  {/if}
</div>

<style>
  .peak-table-container {
    height: 100%;
    display: flex;
    flex-direction: column;
  }

  .peak-table-container.carousel {
    height: 100%;
    max-height: none;
  }

  .table-wrapper {
    flex: 1;
    overflow: auto;
    scrollbar-width: thin;
    border-radius: var(--enforce-pill);
    background: rgba(255, 255, 255, 0.85);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04) inset;
    height: 100%;
    max-height: 100%;
  }

  .peak-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.85rem;
  }

  th {
    position: sticky;
    top: 0;
    background-color: rgba(255, 255, 255, 0.95);
    backdrop-filter: blur(16px);
    font-weight: 500;
    text-align: left;
    padding: 0.85rem 1.25rem;
    user-select: none;
    cursor: pointer;
    color: var(--text-secondary);
    border-bottom: 1px solid var(--surface-stroke);
    transition: background-color var(--transition-smooth);
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    z-index: 10;
  }

  th:first-child {
    border-top-left-radius: var(--enforce-sm);
  }

  th:last-child {
    border-top-right-radius: var(--enforce-sm);
  }

  th:hover {
    background-color: rgba(240, 240, 245, 0.9);
  }

  th.active {
    color: var(--accent);
  }

  .sort-indicator {
    margin-left: 0.25rem;
    opacity: 0.7;
  }

  tr:hover {
    background-color: rgba(0, 0, 0, 0.03);
  }

  td {
    padding: 0.75rem 1.25rem;
    border-bottom: 1px solid rgba(0, 0, 0, 0.05);
  }

  .intensity-cell {
    position: relative;
    width: 18%;
  }

  .intensity-bar {
    position: relative;
    display: block;
    width: 100%;
    z-index: 1;
  }

  .intensity-bar::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    height: 100%;
    width: var(--percentage);
    background: linear-gradient(90deg, var(--accent-soft) 0%, transparent 100%);
    z-index: -1;
    border-radius: 8px;
  }

  .empty-state {
    display: flex;
    align-items: center;
    justify-content: center;
    height: 100%;
    color: var(--text-secondary);
  }

  .empty-content {
    display: flex;
    flex-direction: column;
    align-items: center;
    opacity: 0.6;
  }

  .empty-content :global(.empty-icon) {
    margin-bottom: 1rem;
    opacity: 0.5;
  }

  .table-wrapper::-webkit-scrollbar {
    width: 4px;
  }

  .table-wrapper::-webkit-scrollbar-thumb {
    background-color: rgba(0, 0, 0, 0.15);
    border-radius: 2px;
  }

  .table-wrapper::-webkit-scrollbar-track {
    background-color: rgba(255, 255, 255, 0.3);
  }
</style>
