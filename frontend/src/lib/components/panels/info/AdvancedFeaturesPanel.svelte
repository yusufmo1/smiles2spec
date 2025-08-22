<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
  import { fade, fly } from 'svelte/transition';
  import { onMount } from 'svelte';

  let hoveredFeature: string | null = null;

  const features = [
    {
      id: 'bulk',
      icon: 'UploadIcon',
      title: 'Bulk Processing',
      description: 'Process hundreds of molecules at once',
      details: 'Upload CSV files with SMILES structures',
      stats: '1000+ molecules/batch',
      color: 'var(--accent)',
    },
    {
      id: 'chat',
      icon: 'ChatWithSpectrumIcon',
      title: 'AI Chat Analysis',
      description: 'Interactive spectrum interpretation',
      details: 'Powered by advanced language models',
      stats: 'Real-time insights',
      color: 'var(--info)',
    },
    {
      id: 'export',
      icon: 'ExportIcon',
      title: 'Multiple Formats',
      description: 'Export in various standard formats',
      details: 'JSON, CSV, MSP, and more',
      stats: '5+ formats',
      color: 'var(--success)',
    },
    {
      id: 'generate',
      icon: 'WandIcon',
      title: 'SMILES Generation',
      description: 'AI-powered molecular generation',
      details: 'Create random valid molecules',
      stats: 'Unlimited creativity',
      color: 'var(--warning)',
    },
  ];

  const additionalFeatures = [
    { icon: 'FlashIcon', label: 'Real-time predictions' },
    { icon: 'GridIcon', label: 'Batch visualization' },
  ];
</script>

<div class="advanced-features">
  <div class="header-section">
    <DynamicIcon name="WandIcon" size={40} color="var(--accent)" />
    <h3>Advanced Features</h3>
    <p class="subtitle">Powerful tools for mass spectrometry analysis</p>
  </div>

  <div class="features-grid">
    {#each features as feature, i}
      <div
        class="feature-card"
        class:hovered={hoveredFeature === feature.id}
        role="button"
        tabindex="0"
        on:mouseenter={() => (hoveredFeature = feature.id)}
        on:mouseleave={() => (hoveredFeature = null)}
        on:focus={() => (hoveredFeature = feature.id)}
        on:blur={() => (hoveredFeature = null)}
        on:keydown={(e) =>
          e.key === 'Enter' && (hoveredFeature = hoveredFeature === feature.id ? null : feature.id)}
        in:fly={{ y: 20, delay: i * 100, duration: 500 }}
      >
        <div class="card-header">
          <div class="icon-wrapper" style="--feature-color: {feature.color}">
            <DynamicIcon name={feature.icon as any} size={32} color={feature.color} />
          </div>
          <div class="feature-stats">{feature.stats}</div>
        </div>

        <div class="card-content">
          <h4>{feature.title}</h4>
          <p class="description">{feature.description}</p>
          {#if hoveredFeature === feature.id}
            <p class="details" in:fade={{ duration: 200 }}>
              {feature.details}
            </p>
          {/if}
        </div>
      </div>
    {/each}
  </div>

  <div class="additional-features" in:fade={{ delay: 500 }}>
    <h4>Also includes:</h4>
    <div class="feature-pills">
      {#each additionalFeatures as feature}
        <div class="feature-pill">
          <DynamicIcon name={feature.icon as any} size={16} />
          <span>{feature.label}</span>
        </div>
      {/each}
    </div>
  </div>
</div>

<style>
  .advanced-features {
    display: flex;
    flex-direction: column;
    gap: 2rem;
    padding: 0.5rem;
  }

  .header-section {
    text-align: center;
  }

  h3 {
    margin: 0.5rem 0;
    font-size: 1.5rem;
    color: var(--text-primary);
    font-weight: 600;
  }

  .subtitle {
    margin: 0;
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .features-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
    gap: 1.5rem;
  }

  .feature-card {
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 12px;
    padding: 1.5rem;
    display: flex;
    flex-direction: column;
    gap: 1rem;
    transition: all 0.3s ease;
    cursor: pointer;
    position: relative;
    overflow: hidden;
  }

  .feature-card::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 3px;
    background: var(--feature-color);
    transform: scaleX(0);
    transition: transform 0.3s ease;
  }

  .feature-card:hover {
    transform: translateY(-4px);
    border-color: var(--feature-color);
    box-shadow: 0 8px 24px rgba(var(--accent-rgb), 0.15);
  }

  .feature-card:hover::before {
    transform: scaleX(1);
  }

  .card-header {
    display: flex;
    justify-content: space-between;
    align-items: flex-start;
  }

  .icon-wrapper {
    width: 60px;
    height: 60px;
    background: linear-gradient(135deg, rgba(var(--accent-rgb), 0.1), transparent);
    border: 1px solid var(--feature-color);
    border-radius: 12px;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: all 0.3s ease;
  }

  .feature-card:hover .icon-wrapper {
    transform: scale(1.1);
    background: linear-gradient(135deg, rgba(var(--accent-rgb), 0.2), transparent);
  }

  .feature-stats {
    font-size: 0.85rem;
    color: var(--text-secondary);
    font-weight: 500;
    padding: 0.25rem 0.75rem;
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 20px;
  }

  .card-content {
    flex: 1;
  }

  .card-content h4 {
    margin: 0 0 0.5rem;
    font-size: 1.1rem;
    color: var(--text-primary);
  }

  .description {
    margin: 0;
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
  }

  .details {
    margin: 0.5rem 0 0;
    font-size: 0.85rem;
    color: var(--text-secondary);
    opacity: 0.8;
  }

  .additional-features {
    text-align: center;
  }

  .additional-features h4 {
    margin: 0 0 1rem;
    font-size: 1.1rem;
    color: var(--text-primary);
  }

  .feature-pills {
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    justify-content: center;
  }

  .feature-pill {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.5rem 1rem;
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 20px;
    font-size: 0.85rem;
    color: var(--text-secondary);
  }

  @media (max-width: 768px) {
    .features-grid {
      grid-template-columns: 1fr;
    }
  }
</style>
