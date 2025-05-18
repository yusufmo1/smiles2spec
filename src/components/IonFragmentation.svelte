<script>
  import { onMount, afterUpdate } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { defaultTheme, defaultConfig, createLightModePlot } from '../services/plotlyTheme.js';
  import { focusedPanel } from '../stores.js';
  
  export let data = null;
  
  let plotElement;
  let topFragments = [];
  let isInCarousel = false;
  
  // Check if component is in carousel when panel is focused
  $: isInCarousel = $focusedPanel !== null;
  
  $: if (data) {
    processData();
  }
  
  $: if (plotElement && topFragments.length > 0) {
    renderPlot();
  }
  
  // Force re-render when this component is shown in the carousel
  $: if ($focusedPanel && plotElement && topFragments.length > 0) {
    // Short delay to ensure DOM is ready
    setTimeout(() => renderPlot(), 100);
  }
  
  // Handle window resize events for the plot
  onMount(() => {
    const handleResize = () => {
      if (plotElement && topFragments.length > 0) {
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
  
  function processData() {
    // Extract the top 10 fragments by intensity
    const fragmentData = [];
    for (let i = 0; i < data.x.length; i++) {
      fragmentData.push({
        mz: data.x[i],
        intensity: data.y[i]
      });
    }
    
    topFragments = fragmentData
      .filter(f => f.intensity > 0.01) // Filter out noise
      .sort((a, b) => b.intensity - a.intensity) // Sort by intensity desc
      .slice(0, 10); // Take top 10
  }
  
  function renderPlot() {
    const xValues = topFragments.map(f => f.mz);
    const yValues = topFragments.map(f => f.intensity);
    
    // Calculate possible fragment formulas based on M/Z
    const annotations = topFragments.map((fragment, index) => {
      const mz = fragment.mz;
      let formula = "";
      
      // Simple approximation for fragment formulas (this would be more complex in a real app)
      if (mz > 100) {
        formula = `C${Math.floor(mz/12)}H${Math.floor(mz/2)}`;
      } else if (mz > 50) {
        formula = `C${Math.floor(mz/12)}H${Math.floor(mz/2)-2}`;
      } else if (mz > 20) {
        formula = `C${Math.floor(mz/12)}H${Math.floor(mz/3)}`;
      } else {
        formula = `H${Math.floor(mz)}`;
      }
      
      return {
        x: fragment.mz,
        y: fragment.intensity,
        text: formula,
        showarrow: false,
        font: {
          size: 11,
          color: 'var(--text-secondary)'
        },
        yshift: 15,
        xshift: 0
      };
    });
    
    const plotData = [{
      x: xValues,
      y: yValues,
      type: 'bar',
      marker: {
        color: 'var(--accent)',
        line: {
          width: 1,
          color: 'rgba(0,0,0,0.05)'
        },
        opacity: 0.9
      },
      text: topFragments.map(f => `${Math.round(f.intensity * 100)}%`),
      textposition: 'auto',
      hovertemplate: 'M/Z: %{x}<br>Intensity: %{y:.4f}<extra></extra>'
    }];
    
    const layout = {
      title: {
        text: 'Top Fragment Ions',
        font: {
          size: 0 // Hide title as it's in the panel header
        }
      },
      xaxis: {
        title: {
          text: 'M/Z',
          font: {
            size: 12,
            color: 'var(--text-secondary)'
          }
        }
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
      annotations: annotations,
      bargap: 0.3,
      autosize: true,
      responsive: true,
      hovermode: 'closest' // Consistent hover behavior
    };
    
    createLightModePlot(plotElement, plotData, layout);
  }
</script>

<div class="plot-wrapper">
  {#if topFragments.length > 0}
    <div bind:this={plotElement} class="plot-container"></div>
  {:else}
    <div class="placeholder">
      <div class="placeholder-content">
        <span class="placeholder-icon">⚛️</span>
        <p>Awaiting fragmentation data</p>
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
    overflow: hidden;
    flex: 1;
    position: relative;
  }
  
  .placeholder {
    width: 100%;
    height: 100%;
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