<script lang="ts">
  // About page component with static content
  // Panel system is available in dedicated routes
</script>

<div class="about-page">
  <div class="body-wrapper">
    <div class="about-content">
      <div class="about-section">
        <h3>About MS Predictor</h3>
        <p>
          MS Predictor is an advanced AI-powered tool for predicting high-resolution mass spectra
          from SMILES notation. Our system combines machine learning algorithms with extensive
          spectral databases to provide accurate and reliable predictions.
        </p>

        <h3>Key Features</h3>
        <p>
          • AI-powered spectrum prediction<br />
          • SMILES notation input support<br />
          • High-resolution LC-ESI-qTOF simulation<br />
          • Bulk processing capabilities<br />
          • Interactive visualization tools<br />
          • Export functionality for research use
        </p>

        <h3>Technology Stack</h3>
        <p>
          Built with modern web technologies including SvelteKit for the frontend, Python-based
          machine learning models for predictions, and advanced visualization libraries for
          interactive spectrum display.
        </p>
      </div>
    </div>
  </div>
</div>

<style>
  .about-page {
    width: 100%;
    padding: 0;
    margin-bottom: 2rem;
  }

  .body-wrapper {
    width: 100%;
  }

  .about-content {
    background: rgba(255, 255, 255, 0.6);
    border-radius: 16px;
    padding: 2rem;
    box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
    backdrop-filter: blur(10px);
  }

  .about-section h3 {
    font-size: 1.25rem;
    color: #7879ff;
    margin-bottom: 0.75rem;
    margin-top: 1.5rem;
  }

  .about-section h3:first-child {
    margin-top: 0;
  }

  .about-section p {
    font-size: 1rem;
    line-height: 1.6;
    color: #666;
    margin: 0 0 1rem 0;
  }

  @media (max-width: 768px) {
    .about-content {
      padding: 1.5rem;
    }

    .about-section h3 {
      font-size: 1.1rem;
    }

    .about-section p {
      font-size: 0.95rem;
    }
  }
</style>
