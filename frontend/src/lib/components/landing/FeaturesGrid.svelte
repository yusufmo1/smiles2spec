<script lang="ts">
  import { DynamicIcon, type IconName } from '$lib/components/icons';

  interface Feature {
    icon: IconName;
    title: string;
    description: string;
    color: string;
  }

  const features: Feature[] = [
    {
      icon: 'FlashIcon',
      title: 'Instant Predictions',
      description:
        'Get accurate EI mass spectra from SMILES strings in <100ms using optimised Random Forest models with 8000+ molecular descriptors.',
      color: 'var(--accent)',
    },
    {
      icon: 'TargetIcon',
      title: 'High Accuracy',
      description:
        '99.2% prediction accuracy achieved through comprehensive training on diverse molecular structures and validated spectral data.',
      color: 'var(--accent-secondary)',
    },
    {
      icon: 'ScikitLearnIcon',
      title: 'Advanced ML Pipeline',
      description:
        'Foundational Differentiable Random Forest models with RDKit molecular fingerprints, descriptors, and custom feature engineering.',
      color: '#f7931e',
    },
    {
      icon: 'DatabaseIcon',
      title: 'Rich Feature Set',
      description:
        'Over 8000 molecular descriptors including topological, electronic, and physicochemical properties for comprehensive analysis.',
      color: '#4fc3f7',
    },
    {
      icon: 'UserIcon',
      title: 'Easy to Use',
      description:
        'No complex setup required - just paste your SMILES and predict. Supports bulk uploads and multiple export formats.',
      color: '#34c759',
    },
    {
      icon: 'HeartIcon',
      title: 'Open Source',
      description:
        'Free for researchers worldwide, built for the scientific community. Full source code available with comprehensive documentation.',
      color: '#ff453a',
    },
  ];
</script>

<section class="features-section">
  <div class="features-header">
    <div class="title-pill">
      <div class="pill-title">Why Choose SMILES2SPEC?</div>
      <div class="pill-subtitle">
        Preview our foundational spectral simulation model - advanced machine learning meets
        molecular science for the future of computational chemistry research
      </div>
    </div>
  </div>

  <div class="features-grid">
    {#each features as feature}
      <div class="feature-card glass-card">
        <div class="feature-icon" style="color: {feature.color}">
          <DynamicIcon name={feature.icon} size={32} color={feature.color} />
        </div>
        <h3 class="feature-title">{feature.title}</h3>
        <p class="feature-description">{feature.description}</p>
      </div>
    {/each}
  </div>
</section>

<style>
  .features-section {
    padding: 2rem 1rem;
    max-width: 1200px;
    margin: 0 auto;
  }

  .features-header {
    text-align: center;
    margin-bottom: 3rem;
  }

  .title-pill {
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
    margin: 0 0 2rem 0;
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

  .title-pill::before {
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

  .title-pill:hover {
    transform: translateY(-2px);
    box-shadow:
      0 12px 40px rgba(120, 121, 255, 0.15),
      inset 0 1px 0 rgba(255, 255, 255, 0.15);
  }

  .title-pill:hover::before {
    opacity: 1;
  }

  .pill-title {
    font-size: clamp(1.5rem, 4vw, 2rem);
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

  .features-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
    gap: 1.5rem;
    margin-top: 2rem;
  }

  .feature-card {
    padding: 2rem;
    border-radius: 20px;
    background: rgba(255, 255, 255, 0.08);
    backdrop-filter: blur(20px);
    -webkit-backdrop-filter: blur(20px);
    border: 1px solid rgba(255, 255, 255, 0.1);
    transition: all var(--transition-smooth);
    text-align: center;
    position: relative;
    overflow: hidden;
  }

  .feature-card::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: linear-gradient(
      135deg,
      rgba(255, 255, 255, 0.1) 0%,
      rgba(255, 255, 255, 0.05) 100%
    );
    border-radius: inherit;
    opacity: 0;
    transition: opacity var(--transition-smooth);
  }

  .feature-card:hover {
    transform: translateY(-4px);
    box-shadow: 0 20px 40px rgba(0, 0, 0, 0.1);
  }

  .feature-card:hover::before {
    opacity: 1;
  }

  .feature-icon {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 60px;
    height: 60px;
    margin: 0 auto 1.5rem;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 16px;
    position: relative;
    z-index: 1;
  }

  .feature-title {
    font-size: 1.25rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 1rem 0;
    position: relative;
    z-index: 1;
  }

  .feature-description {
    font-size: 0.95rem;
    color: var(--text-secondary);
    line-height: 1.6;
    margin: 0;
    position: relative;
    z-index: 1;
  }

  /* Responsive adjustments */
  @media (max-width: 768px) {
    .features-section {
      padding: 2rem 1rem;
    }

    .title-pill {
      padding: 1.5rem 2rem;
      max-width: 90%;
      margin: 0 auto 2rem auto;
    }

    .features-grid {
      grid-template-columns: 1fr;
      gap: 1rem;
    }

    .feature-card {
      padding: 1.5rem;
    }
  }

  @media (max-width: 480px) {
    .title-pill {
      padding: 1.25rem 1.5rem;
      border-radius: 24px;
    }

    .pill-title {
      margin-bottom: 0.5rem;
    }

    .feature-card {
      padding: 1.25rem;
    }

    .feature-icon {
      width: 50px;
      height: 50px;
      margin-bottom: 1rem;
    }
  }
</style>
