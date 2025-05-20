<script>
  // Component for the "How It Works" page
</script>

<div class="how-it-works">
  <div class="glass-card info-card">
    <h2>How MS Predictor Works</h2>
    
    <div class="info-section">
      <h3>SMILES Input</h3>
      <p>Enter a SMILES (Simplified Molecular Input Line Entry System) string to represent the chemical structure you want to analyze.</p>
    </div>
    
    <div class="info-section">
      <h3>Prediction Engine</h3>
      <p>Our advanced AI model processes the SMILES notation and generates predictions for how the molecule would appear in a high-resolution LC-ESI-qTOF mass spectrum.</p>
    </div>
    
    <div class="info-section">
      <h3>Spectral Analysis</h3>
      <p>The system calculates expected fragmentation patterns, peak intensities, and m/z ratios that would be observed in a real-world mass spectrometry experiment.</p>
    </div>
    
    <div class="info-section">
      <h3>Visualization</h3>
      <p>Results are displayed in multiple interactive formats including spectra plots, fragment ion breakdowns, and data tables for comprehensive analysis.</p>
    </div>
  </div>
</div>

<style>
  .how-it-works {
    width: 100%;
    padding: 0 1.5rem;
    margin-bottom: 2rem;
  }
  
  .info-card {
    padding: 2.5rem;
    width: 100%;
  }
  
  h2 {
    font-size: 1.8rem;
    margin: 0 0 2rem;
    text-align: center;
    background: linear-gradient(135deg, #7879ff 0%, #6d56f7 60%, #7879ff 100%) 0 0/200%;
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    animation: title-shift 15s linear infinite;
  }
  
  .info-section {
    margin-bottom: 1.75rem;
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
  
  /* Responsive adjustments */
  @media (max-width: 768px) {
    .info-card {
      padding: 1.75rem;
    }
    
    h2 {
      font-size: 1.5rem;
      margin-bottom: 1.5rem;
    }
    
    h3 {
      font-size: 1.1rem;
    }
  }
</style> 