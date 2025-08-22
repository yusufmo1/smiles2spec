<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
  import { fade, scale } from 'svelte/transition';
  import { onMount } from 'svelte';

  let selectedFeature = 0;
  let demoData = {
    spectrum: [
      { x: 50, y: 15 },
      { x: 75, y: 30 },
      { x: 109, y: 100 },
      { x: 125, y: 45 },
      { x: 150, y: 20 },
      { x: 194, y: 85 },
    ],
    peaks: [
      { mz: 109.05, intensity: 100, formula: 'C5H5N2O' },
      { mz: 194.08, intensity: 85, formula: 'C8H10N4O2' },
      { mz: 125.06, intensity: 45, formula: 'C6H7N2O' },
    ],
  };

  const features = [
    {
      id: 'spectrum',
      icon: 'SpectrumIcon',
      title: 'Interactive Spectrum',
      description: 'Real-time visualization with zoom & pan',
      capabilities: [
        'Zoom into specific m/z ranges',
        'Hover for peak details',
        'Click to select peaks',
        'Export as PNG/SVG',
      ],
    },
    {
      id: 'table',
      icon: 'DataTableIcon',
      title: 'Peak Analysis',
      description: 'Comprehensive data exploration',
      capabilities: [
        'Sort by m/z or intensity',
        'Filter by threshold',
        'Chemical formula display',
        'CSV/JSON export',
      ],
    },
    {
      id: 'structure',
      icon: 'MicroscopeIcon',
      title: 'Structure Viewer',
      description: '2D molecular visualization',
      capabilities: [
        'Fragment highlighting',
        'Stereochemistry display',
        'Atom numbering',
        'Export structure',
      ],
    },
  ];

  function selectFeature(index: number) {
    selectedFeature = index;
  }

  // Simulated spectrum drawing
  function drawSpectrum(canvas: HTMLCanvasElement) {
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const width = canvas.width;
    const height = canvas.height;

    // Clear canvas
    ctx.clearRect(0, 0, width, height);

    // Draw axes
    ctx.strokeStyle = 'var(--glass-border)';
    ctx.beginPath();
    ctx.moveTo(20, height - 20);
    ctx.lineTo(width - 20, height - 20);
    ctx.moveTo(20, 20);
    ctx.lineTo(20, height - 20);
    ctx.stroke();

    // Draw peaks
    ctx.strokeStyle = 'var(--accent)';
    ctx.lineWidth = 2;

    demoData.spectrum.forEach((peak) => {
      const x = 20 + (peak.x / 200) * (width - 40);
      const y = height - 20 - (peak.y / 100) * (height - 40);

      ctx.beginPath();
      ctx.moveTo(x, height - 20);
      ctx.lineTo(x, y);
      ctx.stroke();
    });
  }

  onMount(() => {
    const canvas = document.querySelector('.demo-spectrum') as HTMLCanvasElement;
    if (canvas) {
      drawSpectrum(canvas);
    }
  });
</script>

<div class="visualization-features">
  <div class="header-section">
    <DynamicIcon name="GridIcon" size={40} color="var(--accent)" />
    <h3>Visualization Features</h3>
    <p class="subtitle">Interactive tools for spectrum analysis</p>
  </div>

  <div class="features-showcase">
    <div class="feature-tabs">
      {#each features as feature, i}
        <button
          class="feature-tab"
          class:active={selectedFeature === i}
          on:click={() => selectFeature(i)}
          transition:scale={{ duration: 200 }}
        >
          <DynamicIcon name={feature.icon as any} size={24} />
          <span>{feature.title}</span>
        </button>
      {/each}
    </div>

    <div class="feature-content">
      {#key selectedFeature}
        <div class="content-wrapper" in:fade={{ duration: 300 }}>
          <div class="demo-area">
            {#if features[selectedFeature].id === 'spectrum'}
              <canvas class="demo-spectrum" width="320" height="200"></canvas>
            {:else if features[selectedFeature].id === 'table'}
              <div class="demo-table">
                <table>
                  <thead>
                    <tr>
                      <th>m/z</th>
                      <th>Intensity</th>
                      <th>Formula</th>
                    </tr>
                  </thead>
                  <tbody>
                    {#each demoData.peaks as peak}
                      <tr>
                        <td>{peak.mz.toFixed(2)}</td>
                        <td>{peak.intensity}%</td>
                        <td class="formula">{peak.formula}</td>
                      </tr>
                    {/each}
                  </tbody>
                </table>
              </div>
            {:else if features[selectedFeature].id === 'structure'}
              <div class="demo-structure">
                <div class="molecule-placeholder">
                  <DynamicIcon name="MicroscopeIcon" size={80} color="var(--glass-border)" />
                  <p>2D Structure View</p>
                </div>
              </div>
            {/if}
          </div>

          <div class="feature-details">
            <h4>{features[selectedFeature].title}</h4>
            <p class="description">{features[selectedFeature].description}</p>

            <ul class="capabilities">
              {#each features[selectedFeature].capabilities as capability}
                <li>{capability}</li>
              {/each}
            </ul>
          </div>
        </div>
      {/key}
    </div>
  </div>
</div>

<style>
  .visualization-features {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
    padding: 0.5rem;
  }

  .header-section {
    text-align: center;
  }

  h3 {
    margin: 0.5rem 0;
    font-size: 1.5rem;
    color: var(--text-primary);
    font-weight: 600;
  }

  .subtitle {
    margin: 0;
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .features-showcase {
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 12px;
    overflow: hidden;
  }

  .feature-tabs {
    display: flex;
    border-bottom: 1px solid var(--glass-border);
  }

  .feature-tab {
    flex: 1;
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.5rem;
    padding: 1rem;
    background: transparent;
    border: none;
    color: var(--text-secondary);
    cursor: pointer;
    transition: all 0.3s ease;
    position: relative;
  }

  .feature-tab:hover {
    background: var(--glass-bg);
    color: var(--text-primary);
  }

  .feature-tab.active {
    color: var(--accent);
    background: var(--glass-bg);
  }

  .feature-tab.active::after {
    content: '';
    position: absolute;
    bottom: -1px;
    left: 0;
    right: 0;
    height: 2px;
    background: var(--accent);
  }

  .feature-tab span {
    font-size: 0.9rem;
    font-weight: 500;
  }

  .feature-content {
    padding: 1.5rem;
  }

  .content-wrapper {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 2rem;
    align-items: center;
  }

  .demo-area {
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 200px;
  }

  .demo-spectrum {
    border: 1px solid var(--glass-border);
    border-radius: 8px;
    background: var(--surface);
  }

  .demo-table {
    width: 100%;
    max-width: 320px;
  }

  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.85rem;
  }

  th {
    text-align: left;
    padding: 0.5rem;
    border-bottom: 2px solid var(--glass-border);
    color: var(--text-primary);
    font-weight: 600;
  }

  td {
    padding: 0.5rem;
    border-bottom: 1px solid var(--glass-border);
    color: var(--text-secondary);
  }

  .formula {
    font-family: var(--font-mono);
    color: var(--accent);
  }

  .demo-structure {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 320px;
    height: 200px;
    border: 1px solid var(--glass-border);
    border-radius: 8px;
    background: var(--surface);
  }

  .molecule-placeholder {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0.5rem;
    color: var(--text-secondary);
  }

  .molecule-placeholder p {
    margin: 0;
    font-size: 0.9rem;
  }

  .feature-details h4 {
    margin: 0 0 0.5rem;
    font-size: 1.2rem;
    color: var(--text-primary);
  }

  .description {
    margin: 0 0 1rem;
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .capabilities {
    list-style: none;
    padding: 0;
    margin: 0;
  }

  .capabilities li {
    position: relative;
    padding-left: 1.5rem;
    margin-bottom: 0.5rem;
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .capabilities li::before {
    content: '✓';
    position: absolute;
    left: 0;
    color: var(--accent);
  }

  @media (max-width: 768px) {
    .content-wrapper {
      grid-template-columns: 1fr;
      gap: 1.5rem;
    }

    .feature-tabs {
      flex-wrap: wrap;
    }

    .feature-tab span {
      display: none;
    }
  }
</style>
