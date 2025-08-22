<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
  import BentoGrid from './components/BentoGrid.svelte';

  const stats = [
    {
      id: 'accuracy',
      title: 'Accuracy',
      value: '95%+',
      icon: 'TargetIcon',
      color: '#34c759',
      size: 'medium' as const,
    },
    {
      id: 'speed',
      title: 'Speed',
      value: '<2s',
      icon: 'FlashIcon',
      color: '#ff9f0a',
      size: 'small' as const,
    },
    {
      id: 'range',
      title: 'MW Range',
      value: '1000 Da',
      icon: 'ScaleIcon',
      color: '#5856d6',
      size: 'small' as const,
    },
    {
      id: 'spectra',
      title: 'Training Data',
      value: '50K+',
      icon: 'SpectrumIcon',
      color: '#7879ff',
      size: 'medium' as const,
    },
  ];
</script>

<div class="overview-panel">
  <div class="hero-section">
    <div class="hero-icon">
      <div class="icon-pulse"></div>
      <DynamicIcon name="SpectrumIcon" size={48} color="white" />
    </div>
    <h1 class="hero-title">AI-POWERED PREDICTION</h1>
    <p class="hero-subtitle">Transform SMILES into accurate mass spectra in seconds</p>
  </div>

  <div class="stats-section">
    <BentoGrid items={stats} />
  </div>
</div>

<style>
  .overview-panel {
    height: 100%;
    display: flex;
    flex-direction: column;
    gap: 2rem;
  }

  .hero-section {
    text-align: center;
    padding: 2rem 1rem 1rem;
  }

  .hero-icon {
    width: 100px;
    height: 100px;
    background: linear-gradient(135deg, var(--accent), var(--accent-secondary));
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    margin: 0 auto 1.5rem;
    position: relative;
    box-shadow: 0 12px 32px rgba(120, 121, 255, 0.3);
  }

  .icon-pulse {
    position: absolute;
    inset: -20px;
    background: radial-gradient(circle, var(--accent), transparent 70%);
    opacity: 0.3;
    animation: pulse 3s ease-in-out infinite;
  }

  @keyframes pulse {
    0%,
    100% {
      transform: scale(1);
      opacity: 0.3;
    }
    50% {
      transform: scale(1.2);
      opacity: 0.1;
    }
  }

  .hero-title {
    font-size: 1.75rem;
    font-weight: 700;
    margin: 0 0 0.5rem;
    background: linear-gradient(135deg, var(--accent), var(--accent-secondary));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    letter-spacing: 0.05em;
  }

  .hero-subtitle {
    font-size: 1rem;
    color: var(--text-secondary);
    margin: 0;
    max-width: 300px;
    margin: 0 auto;
  }

  .stats-section {
    flex: 1;
  }
</style>
