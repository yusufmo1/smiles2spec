<script lang="ts">
  import { goto } from '$app/navigation';
  import { DynamicIcon } from '$lib/components/icons';

  function navigateToSimulation() {
    goto('/spectral-simulation');
  }

  function navigateToHowItWorks() {
    goto('/how-it-works');
  }
</script>

<section class="hero-section">
  <div class="hero-content">
    <!-- Main Title -->
    <h1 class="hero-title prosto-one-regular">
      <span class="text-gradient">SMILES2SPEC</span>
    </h1>

    <!-- Tagline Pill -->
    <div class="tagline-pill">
      <div class="pill-title">AI-Powered Mass Spectrum Prediction</div>
      <div class="pill-subtitle">
        Transform SMILES molecular notation into high-resolution electron ionization mass spectra
        using advanced machine learning. Preview our foundational spectral simulation model with
        RDKit molecular processing and Differentiable Random Forest algorithms trained on 8000+
        molecular descriptors.
      </div>
    </div>

    <!-- Technical Highlights -->
    <div class="hero-highlights">
      <div class="highlight-item">
        <DynamicIcon name="CpuIcon" size={16} color="var(--accent)" />
        <span>Foundational DRF</span>
      </div>
      <div class="highlight-item">
        <DynamicIcon name="DatabaseIcon" size={16} color="var(--accent-secondary)" />
        <span>8000+ RDKit Descriptors</span>
      </div>
      <div class="highlight-item">
        <DynamicIcon name="FlashIcon" size={16} color="#34c759" />
        <span>Sub-second Predictions</span>
      </div>
      <div class="highlight-item">
        <DynamicIcon name="TargetIcon" size={16} color="#ff453a" />
        <span>99.2% Accuracy</span>
      </div>
    </div>

    <!-- CTA Buttons -->
    <div class="hero-actions">
      <button class="btn btn-primary" on:click={navigateToSimulation}>
        <DynamicIcon name="SpectrumIcon" size={20} />
        Try Preview Model
      </button>

      <button class="btn btn-secondary" on:click={navigateToHowItWorks}>
        <DynamicIcon name="BookIcon" size={20} />
        Learn How It Works
      </button>
    </div>

    <!-- Visual Element -->
    <div class="hero-visual">
      <div class="molecule-animation">
        <DynamicIcon name="MoleculeIcon" size={60} color="var(--accent)" />
      </div>
      <div class="spectrum-preview">
        <DynamicIcon name="ChartIcon" size={50} color="var(--accent-secondary)" />
      </div>
    </div>
  </div>
</section>

<style>
  .hero-section {
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 60vh;
    text-align: center;
    padding: 3rem 1rem;
  }

  .hero-content {
    max-width: 800px;
    width: 100%;
  }

  .hero-title {
    font-size: clamp(2.5rem, 8vw, 4rem);
    font-weight: 400;
    margin: 0 0 1rem 0;
    line-height: 1.1;
  }

  .text-gradient {
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    -webkit-background-clip: text;
    background-clip: text;
    -webkit-text-fill-color: transparent;
    color: transparent;
  }

  .tagline-pill {
    display: inline-flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
    padding: 2rem 3rem;
    background: linear-gradient(
      135deg,
      rgba(120, 121, 255, 0.08) 0%,
      rgba(164, 93, 230, 0.08) 100%
    );
    border: 1px solid rgba(120, 121, 255, 0.2);
    border-radius: 32px;
    margin: 0 auto 2rem auto;
    backdrop-filter: blur(20px);
    -webkit-backdrop-filter: blur(20px);
    transition: all var(--transition-smooth);
    position: relative;
    overflow: hidden;
    max-width: 700px;
    box-shadow:
      0 8px 32px rgba(120, 121, 255, 0.1),
      inset 0 1px 0 rgba(255, 255, 255, 0.1);
  }

  .tagline-pill::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: linear-gradient(
      135deg,
      rgba(120, 121, 255, 0.05) 0%,
      rgba(164, 93, 230, 0.05) 100%
    );
    border-radius: inherit;
    opacity: 0;
    transition: opacity var(--transition-smooth);
  }

  .tagline-pill:hover {
    transform: translateY(-2px);
    box-shadow:
      0 12px 40px rgba(120, 121, 255, 0.15),
      inset 0 1px 0 rgba(255, 255, 255, 0.15);
  }

  .tagline-pill:hover::before {
    opacity: 1;
  }

  .pill-title {
    font-size: clamp(1.25rem, 4vw, 1.75rem);
    font-weight: 700;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    -webkit-background-clip: text;
    background-clip: text;
    -webkit-text-fill-color: transparent;
    color: transparent;
    margin-bottom: 0.75rem;
    line-height: 1.2;
    position: relative;
    z-index: 1;
  }

  .pill-subtitle {
    font-size: clamp(0.95rem, 2.2vw, 1.1rem);
    color: var(--text-secondary);
    line-height: 1.5;
    font-weight: 400;
    position: relative;
    z-index: 1;
  }

  .hero-highlights {
    display: flex;
    justify-content: center;
    gap: 1.5rem;
    flex-wrap: wrap;
    margin: 0 0 2.5rem 0;
  }

  .highlight-item {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.5rem 1rem;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 20px;
    backdrop-filter: blur(10px);
    -webkit-backdrop-filter: blur(10px);
    border: 1px solid rgba(255, 255, 255, 0.1);
    font-size: 0.875rem;
    font-weight: 500;
    color: var(--text-primary);
    transition: all var(--transition-smooth);
  }

  .highlight-item:hover {
    background: rgba(255, 255, 255, 0.12);
    transform: translateY(-1px);
  }

  .hero-actions {
    display: flex;
    gap: 1rem;
    justify-content: center;
    flex-wrap: wrap;
    margin-bottom: 3rem;
  }

  .btn {
    display: inline-flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.875rem 2rem;
    border-radius: 50px;
    font-size: 1rem;
    font-weight: 500;
    text-decoration: none;
    border: none;
    cursor: pointer;
    transition: all var(--transition-smooth);
    backdrop-filter: blur(10px);
    -webkit-backdrop-filter: blur(10px);
  }

  .btn-primary {
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    color: white;
    box-shadow: 0 8px 24px rgba(120, 121, 255, 0.3);
  }

  .btn-primary:hover {
    transform: translateY(-2px);
    box-shadow: 0 12px 32px rgba(120, 121, 255, 0.4);
  }

  .btn-secondary {
    background: rgba(255, 255, 255, 0.1);
    color: var(--text-primary);
    border: 1px solid rgba(255, 255, 255, 0.2);
  }

  .btn-secondary:hover {
    background: rgba(255, 255, 255, 0.15);
    transform: translateY(-2px);
  }

  .hero-visual {
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 2rem;
    opacity: 0.7;
  }

  .molecule-animation,
  .spectrum-preview {
    animation: float 3s ease-in-out infinite;
  }

  .spectrum-preview {
    animation-delay: 1.5s;
  }

  @keyframes float {
    0%,
    100% {
      transform: translateY(0px);
    }
    50% {
      transform: translateY(-10px);
    }
  }

  /* Responsive adjustments */
  @media (max-width: 768px) {
    .tagline-pill {
      padding: 1.5rem 2rem;
      max-width: 90%;
      margin: 0 auto 2rem auto;
    }
  }

  @media (max-width: 640px) {
    .hero-section {
      min-height: 50vh;
      padding: 2rem 1rem;
    }

    .tagline-pill {
      padding: 1.25rem 1.5rem;
      border-radius: 24px;
    }

    .pill-title {
      margin-bottom: 0.5rem;
    }

    .hero-actions {
      flex-direction: column;
      align-items: center;
    }

    .btn {
      width: 100%;
      max-width: 280px;
    }

    .hero-visual {
      gap: 1rem;
    }

    .hero-highlights {
      gap: 0.75rem;
    }

    .highlight-item {
      font-size: 0.8rem;
      padding: 0.4rem 0.8rem;
    }
  }
</style>
