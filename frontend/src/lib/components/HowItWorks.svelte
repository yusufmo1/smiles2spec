<script lang="ts">
  // How It Works component with static content
  // Panel system is available in dedicated routes
</script>

<div class="how-it-works">
  <div class="body-wrapper">
    <div class="how-it-works-content">
      <div class="step">
        <h3>1. Input SMILES</h3>
        <p>
          Enter a SMILES (Simplified Molecular Input Line Entry System) string representing your
          chemical compound. SMILES notation is a standard way to represent molecular structures as
          text.
        </p>
      </div>

      <div class="step">
        <h3>2. AI Prediction</h3>
        <p>
          Our advanced machine learning model analyses the molecular structure and predicts the mass
          spectrum. The AI has been trained on thousands of high-resolution LC-ESI-qTOF spectra.
        </p>
      </div>

      <div class="step">
        <h3>3. Visualization</h3>
        <p>
          View the predicted spectrum with interactive visualization tools. Explore peak
          assignments, intensities, and molecular fragmentation patterns in detail.
        </p>
      </div>

      <div class="step">
        <h3>4. Analysis & Export</h3>
        <p>
          Analyse the results with our built-in tools and export data in various formats for further
          research and documentation purposes.
        </p>
      </div>
    </div>
  </div>
</div>

<style>
  .how-it-works {
    width: 100%;
    padding: 0;
    margin-bottom: 2rem;
  }

  .body-wrapper {
    width: 100%;
  }

  .how-it-works-content {
    display: grid;
    gap: 1.5rem;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
  }

  .step {
    background: rgba(255, 255, 255, 0.6);
    border-radius: 16px;
    padding: 1.5rem;
    box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
    backdrop-filter: blur(10px);
    transition:
      transform 0.2s ease,
      box-shadow 0.2s ease;
  }

  .step:hover {
    transform: translateY(-2px);
    box-shadow: 0 12px 40px rgba(0, 0, 0, 0.15);
  }

  .step h3 {
    font-size: 1.25rem;
    color: #7879ff;
    margin-bottom: 0.75rem;
    margin-top: 0;
  }

  .step p {
    font-size: 1rem;
    line-height: 1.6;
    color: #666;
    margin: 0;
  }

  @media (max-width: 768px) {
    .how-it-works-content {
      grid-template-columns: 1fr;
      gap: 1rem;
    }

    .step {
      padding: 1.25rem;
    }

    .step h3 {
      font-size: 1.1rem;
    }

    .step p {
      font-size: 0.95rem;
    }
  }
</style>
