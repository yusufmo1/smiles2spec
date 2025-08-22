<script lang="ts">
  /**
   * SpectrumPlot Component
   *
   * Interactive mass spectrum visualization using Plotly.js.
   * Displays predicted mass spectra as bar charts with m/z values
   * on x-axis and relative intensities on y-axis.
   *
   * Features:
   * - Automatic filtering of low-intensity peaks
   * - Responsive resizing with plot manager
   * - Hover tooltips with precise values
   * - Carousel mode support for compact display
   * - Placeholder state with instructions
   *
   * @component
   */
  import { onMount, afterUpdate, onDestroy } from 'svelte';
  import { browser } from '$app/environment';
  import { createPlot, resizePlot } from '$lib/services/plotlyService';
  import { plotManager } from '$lib/services/plotManager';
  import { PlaceholderState } from '../components';

  /**
   * @prop {any} spectrumData - Spectrum data with x (m/z) and y (intensity) arrays
   * @prop {boolean} isCarousel - Whether component is in carousel mode
   */
  export let spectrumData: any = null;
  export let isCarousel = false;

  let plotElement: HTMLElement;
  let isPlotReady = false;
  let cleanupFunction: (() => void) | null = null;

  // Reactive: Re-render when data changes and plot is ready
  $: if (browser && plotElement && spectrumData && isPlotReady) {
    renderSpectrum();
  }

  // NEW: Reactive resize when isCarousel changes
  $: if (isCarousel && plotElement && spectrumData) {
    setTimeout(() => {
      resizePlot(plotElement);
    }, 100);
  }

  /**
   * Force plot resize - useful when container dimensions change
   * @public
   */
  export function forceResize() {
    if (plotElement) {
      resizePlot(plotElement);
    }
  }

  /**
   * Render mass spectrum visualization
   *
   * Processes spectrum data, filters low-intensity peaks,
   * and creates interactive bar chart visualization.
   *
   * @private
   */
  async function renderSpectrum() {
    if (!browser || !spectrumData || !plotElement) {
      console.warn('Cannot render spectrum: missing requirements');
      return;
    }

    try {
      // Filter out low-intensity peaks for cleaner visualization
      const filteredData = spectrumData.x
        .map((x: number, i: number) => ({ x, y: spectrumData.y[i] }))
        .filter((d: any) => d.y > 0.01)
        .sort((a: any, b: any) => a.x - b.x);

      if (filteredData.length === 0) {
        console.warn('No significant peaks to display');
        return;
      }

      const plotData = [
        {
          x: filteredData.map((d: any) => d.x),
          y: filteredData.map((d: any) => d.y),
          type: 'bar',
          name: 'Mass Spectrum',
          marker: {
            color: '#7879ff',
            line: {
              color: '#6d56f7',
              width: 1,
            },
            opacity: 0.8,
          },
          hovertemplate: 'M/Z: %{x:.2f}<br>Intensity: %{y:.4f}<extra></extra>',
        },
      ];

      const layout = {
        title: {
          text: '',
          font: { size: 0 }, // Hide title as it's in panel header
        },
        xaxis: {
          title: {
            text: 'M/Z',
            font: { size: 13, color: '#8e8e93' },
          },
          tickfont: { size: 11, color: '#8e8e93' },
          gridcolor: 'rgba(142, 142, 147, 0.1)',
          zerolinecolor: 'rgba(142, 142, 147, 0.2)',
        },
        yaxis: {
          title: {
            text: 'Relative Intensity',
            font: { size: 13, color: '#8e8e93' },
          },
          tickfont: { size: 11, color: '#8e8e93' },
          gridcolor: 'rgba(142, 142, 147, 0.1)',
          zerolinecolor: 'rgba(142, 142, 147, 0.2)',
        },
        showlegend: false,
        hovermode: 'closest',
        bargap: 0.2,
      };

      const config = {
        responsive: true,
        displayModeBar: false,
        staticPlot: false,
      };

      await createPlot(plotElement, plotData, layout, config);
      console.log('Spectrum plot rendered successfully');
    } catch (error) {
      console.error('Failed to render spectrum:', error);
    }
  }

  onMount(() => {
    if (browser) {
      console.log('SpectrumPlot mounted');
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
  {#if spectrumData}
    <div bind:this={plotElement} class="plot-container"></div>
  {:else}
    <PlaceholderState
      icon="SpectrumIcon"
      message="Awaiting spectral data"
      hint="Enter a SMILES string and click 'Predict' to generate mass spectrum"
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
