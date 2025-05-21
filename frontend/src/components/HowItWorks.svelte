<script>
  import Panel from './Panel.svelte';
</script>

<div class="how-it-works">
  <div class="body-wrapper">
    <!-- First row with two panels -->
    <div class="row">
      <div class="col-half">
        <Panel title="SPECTRAL SIMULATION OVERVIEW">
          <div class="overview-content">
            <p>Our advanced AI-powered tool processes SMILES notations to generate high-resolution mass spectra predictions. Learn how each component works below.</p>
          </div>
        </Panel>
      </div>
      <div class="col-half">
        <Panel title="INPUT PROCESSING">
          <div class="info-section">
            <h3>SMILES Input</h3>
            <p>The system accepts SMILES (Simplified Molecular Input Line Entry System) notation, a standardized way to represent chemical structures with ASCII strings.</p>
            
            <div class="example-box">
              <span class="label">Example SMILES</span>
              <code>CC(=O)OC1=CC=CC=C1C(=O)O</code>
              <span class="caption">Represents acetylsalicylic acid (aspirin)</span>
            </div>
          </div>
        </Panel>
      </div>
    </div>
    
    <!-- Second row with two panels -->
    <div class="row">
      <div class="col-half">
        <Panel title="MOLECULAR ANALYSIS">
          <div class="info-section">
            <h3>Structure Recognition</h3>
            <p>Our algorithm parses the SMILES notation to identify molecular patterns, functional groups, and bond types critical for mass spectrometry behavior.</p>
            
            <div class="visualization-placeholder">
              <span class="icon">⚛️</span>
              <p>Visual representation of molecular parsing</p>
            </div>
          </div>
        </Panel>
      </div>
      <div class="col-half">
        <Panel title="FRAGMENTATION PREDICTION">
          <div class="info-section">
            <h3>Mass Spectrum Generation</h3>
            <p>The system calculates the most likely fragmentation patterns and relative intensities based on:</p>
            <ul>
              <li>Bond dissociation energies</li>
              <li>Structural stability of fragments</li>
              <li>Common fragmentation pathways</li>
              <li>Known mass spectrometry behaviors</li>
            </ul>
          </div>
        </Panel>
      </div>
    </div>
    
    <!-- Third row with two panels -->
    <div class="row">
      <div class="col-half">
        <Panel title="VISUALIZATION">
          <div class="info-section">
            <h3>Interactive Outputs</h3>
            <p>Results are displayed in multiple formats:</p>
            <div class="visualization-grid">
              <div class="viz-item">
                <span class="viz-icon">📊</span>
                <span class="viz-label">Mass Spectrum</span>
              </div>
              <div class="viz-item">
                <span class="viz-icon">⚛️</span>
                <span class="viz-label">Structure View</span>
              </div>
              <div class="viz-item">
                <span class="viz-icon">📋</span>
                <span class="viz-label">Peak Data</span>
              </div>
              <div class="viz-item">
                <span class="viz-icon">💬</span>
                <span class="viz-label">Analysis Chat</span>
              </div>
            </div>
          </div>
        </Panel>
      </div>
      <div class="col-half">
        <Panel title="TECHNICAL SPECIFICATIONS">
          <div class="tech-specs">
            <div class="spec-item">
              <h4>Prediction Engine</h4>
              <p>Based on machine learning models trained on experimental data from over 10,000 compounds.</p>
            </div>
            <div class="spec-item">
              <h4>Mass Range</h4>
              <p>Accurately predicts fragments from m/z 50 to 1200 with peak intensity relative to base peak.</p>
            </div>
          </div>
        </Panel>
      </div>
    </div>
    
    <!-- Fourth row with two panels -->
    <div class="row">
      <div class="col-half">
        <Panel title="ACCURACY METRICS">
          <div class="info-section">
            <h3>Performance Analysis</h3>
            <p>Our system achieves industry-leading prediction accuracy:</p>
            <ul>
              <li>Typical m/z prediction accuracy of ±0.1 Da</li>
              <li>Relative intensity correlation of R² > 0.85 compared to experimental data</li>
              <li>Successful prediction for >95% of valid input structures</li>
            </ul>
          </div>
        </Panel>
      </div>
      <div class="col-half">
        <Panel title="MOLECULE COMPATIBILITY">
          <div class="info-section">
            <h3>Supported Structures</h3>
            <p>Works best with small to medium organic molecules:</p>
            <ul>
              <li>Molecular weight &lt; 1000 Da</li>
              <li>Contains C, H, N, O, P, S, F, Cl, Br, and I atoms</li>
              <li>Single, double, and triple bonds</li>
              <li>Common functional groups and ring structures</li>
            </ul>
          </div>
        </Panel>
      </div>
    </div>
  </div>
</div>

<style>
  /* Add the same grid layout CSS as the main page */
  .how-it-works {
    width: 100%;
    padding: 0;
    margin-bottom: 2rem;
  }
  
  .body-wrapper {
    display: flex;
    flex-direction: column;
    gap: 2rem;
  }
  
  .row {
    display: flex;
    flex-wrap: wrap;
    margin: 0;
  }
  
  .col-half {
    width: 50%;
    padding: 0 1.5rem;
    margin-bottom: 0;
  }
  
  /* Additional styling for the new content elements */
  .info-section {
    height: 100%;
    display: flex;
    flex-direction: column;
  }
  
  h3 {
    font-size: 1.25rem;
    color: var(--accent);
    margin-bottom: 0.75rem;
  }
  
  p {
    font-size: 1rem;
    line-height: 1.6;
    color: var(--text-secondary);
    margin: 0;
  }
  
  ul {
    padding-left: 1.25rem;
    margin: 0.5rem 0 0 0;
  }
  
  li {
    color: var(--text-secondary);
    margin-bottom: 0.5rem;
    line-height: 1.5;
  }
  
  .example-box {
    margin-top: 1.5rem;
    padding: 1rem;
    background: rgba(0, 0, 0, 0.04);
    border-radius: 12px;
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
  }
  
  .label {
    font-size: 0.8rem;
    color: var(--text-tertiary);
  }
  
  code {
    font-family: 'SF Mono', monospace;
    background: rgba(255, 255, 255, 0.8);
    padding: 0.5rem;
    border-radius: 8px;
    font-size: 0.9rem;
    color: var(--accent);
  }
  
  .caption {
    font-size: 0.8rem;
    font-style: italic;
    color: var(--text-tertiary);
  }
  
  .visualization-placeholder {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    margin-top: 1.5rem;
    background: rgba(255, 255, 255, 0.6);
    padding: 2rem;
    border-radius: 12px;
  }
  
  .icon {
    font-size: 3rem;
    margin-bottom: 1rem;
    color: var(--text-secondary);
  }
  
  .visualization-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1rem;
    margin-top: 1.5rem;
  }
  
  .viz-item {
    background: rgba(255, 255, 255, 0.6);
    padding: 1.5rem;
    border-radius: 12px;
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
  }
  
  .viz-icon {
    font-size: 2rem;
    margin-bottom: 0.5rem;
  }
  
  .viz-label {
    font-size: 0.9rem;
    color: var(--text-secondary);
    font-weight: 500;
  }
  
  .tech-specs {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1.5rem;
  }
  
  .spec-item h4 {
    color: var(--accent);
    margin: 0 0 0.75rem 0;
    font-size: 1.1rem;
  }
  
  /* Responsive adjustments to match main page */
  @media (max-width: 1024px) {
    .col-half {
      width: 100%;
      margin-bottom: 2rem;
    }
    
    .row .col-half:last-child {
      margin-bottom: 0;
    }
    
    .row {
      margin-bottom: 0;
    }
    
    .col-half {
      padding: 0 1rem;
    }
    
    .tech-specs {
      grid-template-columns: 1fr;
    }
  }
  
  @media (max-width: 640px) {
    .visualization-grid {
      grid-template-columns: 1fr;
    }
  }
</style> 