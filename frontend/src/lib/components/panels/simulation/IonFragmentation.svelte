<script lang="ts">
  import { onMount, onDestroy } from 'svelte';
  import { browser } from '$app/environment';
  import { createPlot, resizePlot } from '$lib/services/plotlyService';
  import { plotManager } from '$lib/services/plotManager';
  import { PlaceholderState } from '../components';

  export let data: any = null;
  export let isCarousel = false;

  let plotElement: HTMLDivElement;
  let topFragments: Array<{ mz: number; intensity: number }> = [];
  let isPlotReady = false;
  let cleanupFunction: (() => void) | null = null;

  // Process data when it changes
  $: if (data) {
    processFragmentData();
  }

  // Render plot when ready
  $: if (browser && plotElement && topFragments.length > 0 && isPlotReady) {
    renderFragmentPlot();
  }

  // NEW: Reactive resize when isCarousel changes
  $: if (isCarousel && plotElement && topFragments.length > 0) {
    setTimeout(() => {
      resizePlot(plotElement);
    }, 100);
  }

  // NEW: Export function to manually trigger resize
  export function forceResize() {
    if (plotElement) {
      resizePlot(plotElement);
    }
  }

  function processFragmentData() {
    if (!data || !data.x || !data.y) {
      topFragments = [];
      return;
    }

    // Extract and filter fragments
    const fragmentData = [];
    for (let i = 0; i < data.x.length; i++) {
      fragmentData.push({
        mz: data.x[i],
        intensity: data.y[i],
      });
    }

    topFragments = fragmentData
      .filter((f) => f.intensity > 0.01) // Filter out noise
      .sort((a, b) => b.intensity - a.intensity) // Sort by intensity desc
      .slice(0, 10); // Take top 10

    console.log(`🧬 Processed ${topFragments.length} top fragments`);
  }

  async function renderFragmentPlot() {
    if (!browser || !plotElement || topFragments.length === 0) {
      console.warn('Cannot render fragments: missing requirements');
      return;
    }

    try {
      const xValues = topFragments.map((f) => f.mz);
      const yValues = topFragments.map((f) => f.intensity);

      // Generate simple fragment annotations
      const annotations = topFragments.slice(0, 5).map((fragment) => {
        const mz = fragment.mz;
        let formula = '';

        // Simple approximation for fragment formulas
        if (mz > 100) {
          formula = `C${Math.floor(mz / 12)}H${Math.floor(mz / 2)}`;
        } else if (mz > 50) {
          formula = `C${Math.floor(mz / 12)}H${Math.floor(mz / 2) - 2}`;
        } else if (mz > 20) {
          formula = `C${Math.floor(mz / 12)}H${Math.floor(mz / 3)}`;
        } else {
          formula = `H${Math.floor(mz)}`;
        }

        return {
          x: fragment.mz,
          y: fragment.intensity,
          text: formula,
          showarrow: false,
          font: {
            size: 10,
            color: '#8e8e93',
          },
          yshift: 15,
        };
      });

      const plotData = [
        {
          x: xValues,
          y: yValues,
          type: 'bar',
          name: 'Fragment Ions',
          marker: {
            color: '#ff6b6b',
            line: {
              color: '#ff5252',
              width: 1,
            },
            opacity: 0.8,
          },
          text: topFragments.map((f) => `${Math.round(f.intensity * 100)}%`),
          textposition: 'auto',
          hovertemplate: 'M/Z: %{x:.2f}<br>Intensity: %{y:.4f}<br>%{text}<extra></extra>',
        },
      ];

      const layout = {
        title: {
          text: '',
          font: { size: 0 },
        },
        xaxis: {
          title: {
            text: 'M/Z',
            font: { size: 13, color: '#8e8e93' },
          },
          tickfont: { size: 11, color: '#8e8e93' },
          gridcolor: 'rgba(142, 142, 147, 0.1)',
        },
        yaxis: {
          title: {
            text: 'Relative Intensity',
            font: { size: 13, color: '#8e8e93' },
          },
          tickfont: { size: 11, color: '#8e8e93' },
          gridcolor: 'rgba(142, 142, 147, 0.1)',
        },
        annotations: annotations,
        bargap: 0.3,
        showlegend: false,
        hovermode: 'closest',
      };

      const config = {
        responsive: true,
        displayModeBar: false,
      };

      await createPlot(plotElement, plotData, layout, config);
      console.log('Fragment plot rendered successfully');
    } catch (error) {
      console.error('Failed to render fragments:', error);
    }
  }

  onMount(() => {
    if (browser) {
      console.log('🧬 IonFragmentation mounted');
      isPlotReady = true;

      // Register with plot manager
      if (plotElement) {
        plotManager.registerPlot(plotElement);
      }

      // Handle window resize
      const handleResize = () => {
        if (plotElement) {
          resizePlot(plotElement);
        }
      };

      window.addEventListener('resize', handleResize);
      cleanupFunction = () => {
        window.removeEventListener('resize', handleResize);
        if (plotElement) {
          plotManager.unregisterPlot(plotElement);
        }
      };
    }
  });

  onDestroy(() => {
    if (cleanupFunction) {
      cleanupFunction();
    }
  });
</script>

<div class="plot-wrapper" class:carousel={isCarousel}>
  {#if topFragments.length > 0}
    <div bind:this={plotElement} class="plot-container"></div>
  {:else}
    <PlaceholderState
      icon="FragmentationIcon"
      message="Awaiting fragmentation data"
      hint="Fragment ion analysis will appear after spectrum prediction"
    />
  {/if}
</div>

<style>
  .plot-wrapper {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    overflow: hidden;
  }

  .plot-wrapper.carousel {
    height: 100%;
    max-height: none;
  }

  .plot-container {
    width: 100%;
    height: 100%;
    min-height: 300px;
  }
</style>
