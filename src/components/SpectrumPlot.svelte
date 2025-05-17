<script>
  import { onMount, afterUpdate } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { defaultTheme, defaultConfig, createLightModePlot } from '../services/plotlyTheme.js';
  
  export let spectrumData = null;
  let plotElement;
  
  $: if (plotElement && spectrumData) {
    renderPlot();
  }
  
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
      height: 340,
      autosize: true
    };
    
    createLightModePlot(plotElement, data, layout);
  }
</script>

<div>
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
  .plot-container {
    width: 100%;
    height: 340px;
    overflow: hidden;
  }
  
  .placeholder {
    width: 100%;
    height: 340px;
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