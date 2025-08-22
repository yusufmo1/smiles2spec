<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;

  const journeySteps = [
    {
      id: 'biochemistry',
      title: 'Biochemistry Foundation',
      description: 'Built a strong foundation in biochemical principles and laboratory techniques.',
      icon: 'BeakerIcon',
      year: '2020-2023',
    },
    {
      id: 'ml-experience',
      title: 'Production ML Experience',
      description:
        'Gained hands-on experience developing and deploying machine learning models in production environments.',
      icon: 'CpuIcon',
      year: '2023-2024',
    },
    {
      id: 'msc',
      title: 'MSc AI in Biosciences',
      description:
        'Advanced study combining artificial intelligence with biological and chemical applications at Queen Mary University of London.',
      icon: 'GraduationCapIcon',
      year: '2024-2025',
    },
    {
      id: 'smiles2spec',
      title: 'SMILES2SPEC Project',
      description:
        'Developing an open-source platform for predicting mass spectra from molecular structures using machine learning as part of MSc research.',
      icon: 'MoleculeIcon',
      year: '2024-Present',
    },
  ];
</script>

<div class="journey-panel" class:compact>
  <div class="timeline-container">
    <div class="timeline-grid">
      {#each journeySteps as step, i}
        <div class="timeline-step" style="--delay: {i * 0.15}s">
          <div class="step-icon">
            <DynamicIcon name={step.icon as any} size={compact ? 20 : 24} />
          </div>
          <div class="step-content">
            <div class="step-year">{step.year}</div>
            <h4 class="step-title">{step.title}</h4>
            <p class="step-description">{step.description}</p>
          </div>
          {#if i < journeySteps.length - 1}
            <div class="step-connector" style="--connector-delay: {(i + 1) * 0.15 + 0.4}s">
              <DynamicIcon name="ArrowRightIcon" size={16} />
            </div>
          {/if}
        </div>
      {/each}
    </div>
  </div>
</div>

<style>
  .journey-panel {
    padding: 2.5rem 2rem;
    background: var(--surface-glass);
    border-radius: var(--radius-pill);
    border: 1px solid rgba(255, 255, 255, 0.2);
  }

  .journey-panel.compact {
    padding: 2rem 1.5rem;
  }

  .timeline-container {
    position: relative;
    width: 100%;
    max-width: 900px;
    margin: 0 auto;
  }

  .timeline-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 1.5rem;
    position: relative;
  }

  .timeline-step {
    display: flex;
    flex-direction: column;
    align-items: center;
    opacity: 0;
    transform: translateY(20px);
    animation: fadeInUp 0.6s ease-out forwards;
    animation-delay: var(--delay);
    position: relative;
  }

  .step-icon {
    width: 44px;
    height: 44px;
    background: var(--surface-glass);
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    border: 2px solid var(--accent);
    color: var(--accent);
    z-index: 2;
    box-shadow: 0 4px 12px rgba(120, 121, 255, 0.25);
    position: relative;
    margin-bottom: 1.5rem;
    transition: all 0.3s ease;
  }

  .step-icon:hover {
    transform: scale(1.1);
    box-shadow: 0 6px 16px rgba(120, 121, 255, 0.35);
  }

  .step-content {
    text-align: center;
    width: 100%;
  }

  .step-year {
    background: var(--accent-soft);
    color: var(--accent);
    font-size: 0.75rem;
    font-weight: 600;
    padding: 0.25rem 0.75rem;
    border-radius: var(--radius-pill);
    margin-bottom: 0.75rem;
    display: inline-block;
  }

  .step-title {
    font-size: 0.95rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 0.5rem 0;
    line-height: 1.3;
  }

  .step-description {
    font-size: 0.8rem;
    line-height: 1.4;
    color: var(--text-secondary);
    margin: 0;
  }

  .step-connector {
    position: absolute;
    top: 22px;
    right: -24px;
    color: var(--accent);
    opacity: 0;
    animation: fadeInUp 0.5s ease-out forwards;
    animation-delay: var(--connector-delay);
    z-index: 3;
  }

  @keyframes fadeInUp {
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }

  /* Mobile responsive */
  @media (max-width: 768px) {
    .journey-panel {
      padding: 2rem 1rem;
    }

    .timeline-grid {
      grid-template-columns: 1fr;
      gap: 2rem;
      max-width: 300px;
      margin: 0 auto;
    }

    .timeline-step {
      flex-direction: row;
      text-align: left;
      align-items: flex-start;
      gap: 1rem;
    }

    .step-icon {
      margin-bottom: 0;
      flex-shrink: 0;
    }

    .step-content {
      text-align: left;
    }

    .step-connector {
      display: none;
    }

    .step-title {
      font-size: 1rem;
    }

    .step-description {
      font-size: 0.85rem;
    }
  }

  @media (max-width: 480px) {
    .step-title {
      font-size: 0.9rem;
    }

    .step-description {
      font-size: 0.8rem;
    }

    .step-icon {
      width: 36px;
      height: 36px;
    }
  }
</style>
