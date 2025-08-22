<script lang="ts">
  import { DynamicIcon, type IconName } from '$lib/components/icons';

  interface TechnicalHighlight {
    icon: IconName;
    title: string;
    value: string;
    description: string;
    color: string;
  }

  const highlights: TechnicalHighlight[] = [
    {
      icon: 'RDKitIcon',
      title: 'RDKit Integration',
      value: '8000+',
      description: 'Molecular descriptors and fingerprints',
      color: 'var(--accent)',
    },
    {
      icon: 'ScikitLearnIcon',
      title: 'Foundational DRF',
      value: '99.2%',
      description: 'Prediction accuracy on validation set',
      color: 'var(--accent-secondary)',
    },
    {
      icon: 'FlashIcon',
      title: 'Ultra-Fast',
      value: '<100ms',
      description: 'Average prediction response time',
      color: '#34c759',
    },
    {
      icon: 'ExportIcon',
      title: 'Export Ready',
      value: '4',
      description: 'Output formats (CSV, JSON, MSP, Plot)',
      color: '#ff6b35',
    },
  ];
</script>

<section class="tech-capabilities-section">
  <div class="tech-capabilities-header">
    <div class="title-pill">
      <div class="pill-title">Preview Model Highlights</div>
      <div class="pill-subtitle">
        Key technologies and performance metrics demonstrated in our foundational model preview
      </div>
    </div>
  </div>

  <div class="highlights-grid">
    {#each highlights as highlight}
      <div class="highlight-card">
        <div class="highlight-icon" style="color: {highlight.color}">
          <DynamicIcon name={highlight.icon} size={32} color={highlight.color} />
        </div>
        <div class="highlight-value" style="color: {highlight.color}">{highlight.value}</div>
        <h3 class="highlight-title">{highlight.title}</h3>
        <p class="highlight-description">{highlight.description}</p>
      </div>
    {/each}
  </div>
</section>

<style>
  .tech-capabilities-section {
    padding: 2rem 1rem;
    max-width: 1200px;
    margin: 0 auto;
    background: linear-gradient(
      135deg,
      rgba(120, 121, 255, 0.02) 0%,
      rgba(164, 93, 230, 0.02) 100%
    );
    border-radius: 24px;
    position: relative;
    overflow: hidden;
  }

  .tech-capabilities-section::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background:
      radial-gradient(circle at 20% 80%, rgba(120, 121, 255, 0.05) 0%, transparent 50%),
      radial-gradient(circle at 80% 20%, rgba(164, 93, 230, 0.05) 0%, transparent 50%);
    pointer-events: none;
  }

  .tech-capabilities-header {
    text-align: center;
    margin-bottom: 4rem;
    position: relative;
    z-index: 1;
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

  .highlights-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
    gap: 2rem;
    position: relative;
    z-index: 1;
  }

  .highlight-card {
    text-align: center;
    padding: 2rem 1.5rem;
    background: rgba(255, 255, 255, 0.06);
    backdrop-filter: blur(20px);
    -webkit-backdrop-filter: blur(20px);
    border: 1px solid rgba(255, 255, 255, 0.1);
    border-radius: 20px;
    transition: all var(--transition-smooth);
    position: relative;
    overflow: hidden;
  }

  .highlight-card::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: linear-gradient(
      135deg,
      rgba(255, 255, 255, 0.08) 0%,
      rgba(255, 255, 255, 0.04) 100%
    );
    border-radius: inherit;
    opacity: 0;
    transition: opacity var(--transition-smooth);
  }

  .highlight-card:hover {
    transform: translateY(-4px);
    box-shadow: 0 20px 40px rgba(0, 0, 0, 0.1);
  }

  .highlight-card:hover::before {
    opacity: 1;
  }

  .highlight-icon {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 64px;
    height: 64px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 16px;
    margin: 0 auto 1rem;
    position: relative;
    z-index: 1;
  }

  .highlight-value {
    font-size: 2.25rem;
    font-weight: 700;
    margin-bottom: 0.5rem;
    line-height: 1;
    position: relative;
    z-index: 1;
  }

  .highlight-title {
    font-size: 1.125rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 0.75rem 0;
    position: relative;
    z-index: 1;
  }

  .highlight-description {
    font-size: 0.9rem;
    color: var(--text-secondary);
    margin: 0;
    line-height: 1.5;
    position: relative;
    z-index: 1;
  }

  /* Responsive adjustments */
  @media (max-width: 768px) {
    .tech-capabilities-section {
      padding: 2rem 1rem;
    }

    .title-pill {
      padding: 1.5rem 2rem;
      max-width: 90%;
      margin: 0 auto 2rem auto;
    }

    .highlights-grid {
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 1.5rem;
    }

    .highlight-card {
      padding: 1.5rem 1rem;
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

    .highlights-grid {
      grid-template-columns: 1fr;
    }

    .highlight-icon {
      width: 56px;
      height: 56px;
    }

    .highlight-value {
      font-size: 2rem;
    }
  }
</style>
