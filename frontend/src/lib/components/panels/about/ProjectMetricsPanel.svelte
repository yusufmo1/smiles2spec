<script lang="ts">
  import { onMount } from 'svelte';
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;

  const metrics = [
    {
      id: 'spectra',
      label: 'Spectra in Training Data',
      value: '50,000+',
      icon: 'SpectrumIcon',
      color: '#7879ff',
      bgGradient: 'linear-gradient(135deg, #7879ff15, #7879ff25)',
      shadowColor: '120, 121, 255',
    },
    {
      id: 'accuracy',
      label: 'Prediction Accuracy',
      value: '95%+',
      icon: 'ScaleIcon',
      color: '#48bb78',
      bgGradient: 'linear-gradient(135deg, #48bb7815, #48bb7825)',
      shadowColor: '72, 187, 120',
    },
    {
      id: 'speed',
      label: 'Processing Time',
      value: '>100ms',
      icon: 'FlashIcon',
      color: '#ff8a5c',
      bgGradient: 'linear-gradient(135deg, #ff8a5c15, #ff8a5c25)',
      shadowColor: '255, 138, 92',
    },
    {
      id: 'open-source',
      label: 'Open Source',
      value: '100%',
      icon: 'StarIcon',
      color: '#5ccdff',
      bgGradient: 'linear-gradient(135deg, #5ccdff15, #5ccdff25)',
      shadowColor: '92, 205, 255',
    },
  ];

  let elementsVisible = false;

  function startAnimation() {
    elementsVisible = true;
  }

  onMount(() => {
    // Use Intersection Observer to trigger animations when elements are in view
    const observer = new IntersectionObserver(
      (entries) => {
        if (entries[0].isIntersecting) {
          startAnimation();
          observer.disconnect();
        }
      },
      { threshold: 0.2 }
    );

    const container = document.querySelector('.metrics-grid');
    if (container) {
      observer.observe(container);
    }

    return () => {
      observer.disconnect();
    };
  });
</script>

<div class="metrics-panel" class:compact>
  <div class="metrics-grid" class:compact class:animated={elementsVisible}>
    {#each metrics as metric}
      <div
        class="metric-card"
        style="--accent-color: {metric.color}; --bg-gradient: {metric.bgGradient}; --shadow-color: {metric.shadowColor};"
      >
        <div class="metric-icon">
          <DynamicIcon name={metric.icon as any} size={compact ? 30 : 40} />
        </div>

        <div class="metric-value">{metric.value}</div>
        <div class="metric-label">{metric.label}</div>

        <div class="glow-effect"></div>
      </div>
    {/each}
  </div>
</div>

<style>
  .metrics-panel {
    display: flex;
    justify-content: center;
    align-items: center;
    padding: 3rem 2rem;
  }

  .metrics-panel.compact {
    padding: 2rem 1rem;
  }

  .metrics-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 2rem;
    width: 100%;
    max-width: 1200px;
  }

  .metrics-grid.compact {
    gap: 1.5rem;
  }

  .metric-card {
    position: relative;
    background: var(--bg-gradient);
    border-radius: var(--radius-pill);
    padding: 1.5rem;
    border: 1px solid rgba(255, 255, 255, 0.3);
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
    overflow: hidden;
    backdrop-filter: blur(10px);
    transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
    box-shadow:
      0 4px 20px rgba(var(--shadow-color), 0.1),
      inset 0 1px 0 rgba(255, 255, 255, 0.2);
  }

  .metric-card:hover {
    transform: translateY(-8px) scale(1.02);
    box-shadow:
      0 20px 40px rgba(var(--shadow-color), 0.2),
      0 0 0 1px rgba(var(--shadow-color), 0.3),
      inset 0 1px 0 rgba(255, 255, 255, 0.3);
    background: linear-gradient(
      135deg,
      rgba(var(--shadow-color), 0.05),
      rgba(var(--shadow-color), 0.1)
    );
  }

  .metric-icon {
    color: var(--accent-color);
    margin-bottom: 1rem;
    opacity: 0;
    transform: translateY(-20px) rotate(-10deg);
    transition: all 0.6s cubic-bezier(0.4, 0, 0.2, 1);
    filter: drop-shadow(0 4px 8px rgba(var(--shadow-color), 0.3));
  }

  .metrics-grid.animated .metric-icon {
    opacity: 1;
    transform: translateY(0) rotate(0deg);
  }

  .metric-card:hover .metric-icon {
    transform: translateY(-5px) scale(1.1);
    filter: drop-shadow(0 6px 12px rgba(var(--shadow-color), 0.4));
  }

  .metric-value {
    font-size: 2.5rem;
    font-weight: 800;
    color: var(--accent-color);
    margin-bottom: 0.25rem;
    opacity: 0;
    transform: scale(0.6);
    transition: all 0.7s cubic-bezier(0.4, 0, 0.2, 1) 0.2s;
    text-shadow:
      0 2px 4px rgba(0, 0, 0, 0.1),
      0 1px 2px rgba(var(--shadow-color), 0.3);
    letter-spacing: -0.02em;
    filter: contrast(1.1) saturate(1.2);
    position: relative;
  }

  .metric-value::after {
    content: '';
    position: absolute;
    inset: 0;
    background: radial-gradient(ellipse at center, transparent 40%, rgba(255, 255, 255, 0.03) 100%);
    border-radius: var(--radius-md);
    pointer-events: none;
  }

  .metrics-grid.animated .metric-value {
    opacity: 1;
    transform: scale(1);
  }

  .metric-card:hover .metric-value {
    transform: scale(1.05);
    text-shadow:
      0 4px 8px rgba(0, 0, 0, 0.15),
      0 2px 4px rgba(var(--shadow-color), 0.4);
    filter: contrast(1.2) saturate(1.3);
  }

  .metric-label {
    font-size: 0.95rem;
    color: var(--text-secondary);
    opacity: 0;
    transform: translateY(20px);
    transition: all 0.6s cubic-bezier(0.4, 0, 0.2, 1) 0.4s;
    font-weight: 500;
    line-height: 1.3;
  }

  .metrics-grid.animated .metric-label {
    opacity: 1;
    transform: translateY(0);
  }

  .glow-effect {
    position: absolute;
    inset: -2px;
    border-radius: var(--radius-pill);
    background: linear-gradient(45deg, transparent, rgba(var(--shadow-color), 0.1), transparent);
    opacity: 0;
    z-index: -1;
    transition: opacity 0.4s ease;
  }

  .metric-card:hover .glow-effect {
    opacity: 1;
    animation: glow-pulse 2s ease-in-out infinite;
  }

  @keyframes glow-pulse {
    0%,
    100% {
      background: linear-gradient(45deg, transparent, rgba(var(--shadow-color), 0.1), transparent);
    }
    50% {
      background: linear-gradient(
        45deg,
        rgba(var(--shadow-color), 0.05),
        rgba(var(--shadow-color), 0.2),
        rgba(var(--shadow-color), 0.05)
      );
    }
  }

  /* Enhanced text contrast for better visibility */
  .metric-value {
    position: relative;
    overflow: visible;
  }

  .metric-value::before {
    content: '';
    position: absolute;
    inset: -4px;
    background: radial-gradient(
      circle at center,
      rgba(var(--shadow-color), 0.05) 0%,
      transparent 70%
    );
    border-radius: var(--radius-lg);
    z-index: -1;
    opacity: 0.6;
  }

  @media (max-width: 768px) {
    .metrics-grid {
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 1.5rem;
    }

    .metric-card {
      padding: 1.25rem;
    }

    .metric-value {
      font-size: 2rem;
    }

    .metric-label {
      font-size: 0.85rem;
    }

    .metric-icon {
      margin-bottom: 0.75rem;
    }
  }

  @media (max-width: 480px) {
    .metrics-grid {
      grid-template-columns: repeat(2, 1fr);
      gap: 1rem;
    }

    .metric-card {
      padding: 1rem;
    }

    .metric-value {
      font-size: 1.6rem;
    }

    .metric-label {
      font-size: 0.8rem;
    }
  }
</style>
