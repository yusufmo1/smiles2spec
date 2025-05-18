<script>
  import { onMount, afterUpdate } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { defaultTheme, defaultConfig, createLightModePlot } from '../services/plotlyTheme.js';
  import { focusedPanel } from '../stores.js';
  
  export let spectrumData = null;
  let plotElement;
  let isInCarousel = false;
  
  // Check if component is in carousel when panel is focused
  $: isInCarousel = $focusedPanel !== null;
  
  // Force re-render when this component is shown in the carousel
  $: if ($focusedPanel && plotElement && spectrumData) {
    // Short delay to ensure DOM is ready
    setTimeout(() => renderPlot(), 100);
  }
  
  $: if (plotElement && spectrumData) {
    renderPlot();
  }
  
  // Handle window resize events for the plot
  onMount(() => {
    const handleResize = () => {
      if (plotElement && spectrumData) {
        Plotly.relayout(plotElement, {
          autosize: true
        });
      }
    };
    
    window.addEventListener('resize', handleResize);
    
    return () => {
      window.removeEventListener('resize', handleResize);
    };
  });
  
  function renderPlot() {
    // Filter out zero values for cleaner visualization
    const filteredData = [];
    for (let i = 0; i < spectrumData.x.length; i++) {
      if (spectrumData.y[i] > 0.01) { // Only show peaks above a small threshold
        filteredData.push({
          x: spectrumData.x[i],
          y: spectrumData.y[i]
        });
      }
    }
    
    const data = [{
      x: filteredData.map(d => d.x),
      y: filteredData.map(d => d.y),
      type: 'bar',
      marker: {
        color: 'var(--accent)',
        line: {
          width: 1,
          color: 'rgba(0,0,0,0.05)'
        }
      },
      hovertemplate: 'M/Z: %{x}<br>Intensity: %{y:.4f}<extra></extra>'
    }];
    
    const layout = {
      title: {
        text: 'Mass Spectrum',
        font: {
          size: 0, // Hide title as it's in the panel header
          color: 'var(--text-primary)'
        }
      },
      xaxis: {
        title: {
          text: 'M/Z',
          font: {
            size: 12,
            color: 'var(--text-secondary)'
          }
        },
        range: [0, 250]
      },
      yaxis: {
        title: {
          text: 'Relative Intensity',
          font: {
            size: 12,
            color: 'var(--text-secondary)'
          }
        }
      },
      autosize: true,
      responsive: true,
      hovermode: 'closest' // Consistent hover behavior
    };
    
    createLightModePlot(plotElement, data, layout);
  }
</script>

<div class="plot-wrapper">
  {#if spectrumData}
    <div bind:this={plotElement} class="plot-container"></div>
  {:else}
    <div class="placeholder">
      <div class="placeholder-content">
        <span class="placeholder-icon">📊</span>
        <p>Awaiting spectral data</p>
      </div>
    </div>
  {/if}
</div>

<style>
  .plot-wrapper {
    width: 100%;
    height: 100%;
    display: flex;
    flex-direction: column;
  }

  .plot-container {
    width: 100%;
    height: 100%;
    min-height: 340px;
    overflow: hidden;
    flex: 1;
    position: relative;
  }
  
  .placeholder {
    width: 100%;
    height: 100%;
    min-height: 340px;
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: rgba(255, 255, 255, 0.85);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04) inset;
  }
  
  .placeholder-content {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    opacity: 0.6;
  }
  
  .placeholder-icon {
    font-size: 2.25rem;
    margin-bottom: 1rem;
    color: var(--text-secondary);
  }
  
  .placeholder p {
    font-size: 0.95rem;
    color: var(--text-secondary);
    margin: 0;
  }
</style> 