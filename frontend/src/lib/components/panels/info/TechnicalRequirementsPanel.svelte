<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
  import { fade, scale } from 'svelte/transition';
  import { tweened } from 'svelte/motion';
  import { cubicOut } from 'svelte/easing';
  import { onMount } from 'svelte';

  let activeTab = 'compatibility';

  // Create tweened stores at the top level
  const responseTime = tweened(0, { duration: 1000, easing: cubicOut });
  const concurrentUsers = tweened(0, { duration: 1000, easing: cubicOut });
  const uptime = tweened(0, { duration: 1000, easing: cubicOut });

  const performanceMetrics = [
    { value: responseTime, target: 2, label: 'Response Time', unit: 's' },
    { value: concurrentUsers, target: 100, label: 'Concurrent Users', unit: '+' },
    { value: uptime, target: 99.9, label: 'Uptime', unit: '%' },
  ];

  const compatibilityData = {
    supported: [
      { label: 'Molecular Weight', value: '50-1000 Da', icon: 'ScaleIcon' },
      { label: 'Common Elements', value: 'C, H, N, O, P, S', icon: 'AtomIcon' },
      { label: 'Halogens', value: 'F, Cl, Br, I', icon: 'AtomIcon' },
    ],
    limited: [
      { label: 'Organometallics', value: 'Basic support', icon: 'InfoIcon' },
      { label: 'Large Peptides', value: '>1000 Da limited', icon: 'InfoIcon' },
    ],
    unsupported: [
      { label: 'Polymers', value: 'Not supported', icon: 'CloseIcon' },
      { label: 'Mixtures', value: 'Single compound only', icon: 'CloseIcon' },
    ],
  };

  const performanceSpecs = [
    { metric: 'API Response', value: '>100ms', description: 'Average prediction time' },
    { metric: 'Batch Size', value: '1000', description: 'Max molecules per batch' },
    { metric: 'Rate Limit', value: '60/min', description: 'Requests per minute' },
    { metric: 'Storage', value: '7 days', description: 'Result retention' },
  ];

  onMount(() => {
    responseTime.set(2);
    concurrentUsers.set(100);
    uptime.set(99.9);
  });
</script>

<div class="technical-requirements">
  <div class="header-section">
    <DynamicIcon name="FlashIcon" size={40} color="var(--accent)" />
    <h3>Technical Specifications</h3>
    <p class="subtitle">System requirements and performance metrics</p>
  </div>

  <div class="tabs">
    <button
      class="tab"
      class:active={activeTab === 'compatibility'}
      on:click={() => (activeTab = 'compatibility')}
    >
      <DynamicIcon name="ScaleIcon" size={20} />
      <span>Compatibility</span>
    </button>
    <button
      class="tab"
      class:active={activeTab === 'performance'}
      on:click={() => (activeTab = 'performance')}
    >
      <DynamicIcon name="FlashIcon" size={20} />
      <span>Performance</span>
    </button>
  </div>

  <div class="tab-content">
    {#if activeTab === 'compatibility'}
      <div class="compatibility-section" in:fade={{ duration: 300 }}>
        <div class="compatibility-grid">
          <div class="compat-category">
            <h4>
              <DynamicIcon name="CheckIcon" size={20} color="var(--success)" /> Fully Supported
            </h4>
            <div class="compat-items">
              {#each compatibilityData.supported as item}
                <div class="compat-item supported">
                  <DynamicIcon name={item.icon as any} size={20} />
                  <div>
                    <span class="label">{item.label}</span>
                    <span class="value">{item.value}</span>
                  </div>
                </div>
              {/each}
            </div>
          </div>

          <div class="compat-category">
            <h4>
              <DynamicIcon name="InfoIcon" size={20} color="var(--warning)" /> Limited Support
            </h4>
            <div class="compat-items">
              {#each compatibilityData.limited as item}
                <div class="compat-item limited">
                  <DynamicIcon name={item.icon as any} size={20} />
                  <div>
                    <span class="label">{item.label}</span>
                    <span class="value">{item.value}</span>
                  </div>
                </div>
              {/each}
            </div>
          </div>

          <div class="compat-category">
            <h4><DynamicIcon name="CloseIcon" size={20} color="var(--error)" /> Not Supported</h4>
            <div class="compat-items">
              {#each compatibilityData.unsupported as item}
                <div class="compat-item unsupported">
                  <DynamicIcon name={item.icon as any} size={20} />
                  <div>
                    <span class="label">{item.label}</span>
                    <span class="value">{item.value}</span>
                  </div>
                </div>
              {/each}
            </div>
          </div>
        </div>

        <div class="element-showcase">
          <h5>Supported Elements</h5>
          <div class="element-grid">
            {#each ['C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'] as element}
              <div
                class="element-badge"
                in:scale={{
                  delay: 50 * ['C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'].indexOf(element),
                }}
              >
                {element}
              </div>
            {/each}
          </div>
        </div>
      </div>
    {:else if activeTab === 'performance'}
      <div class="performance-section" in:fade={{ duration: 300 }}>
        <div class="performance-metrics">
          <div class="metric-display" in:fade={{ delay: 0 }}>
            <div class="metric-circle">
              <svg width="100" height="100" viewBox="0 0 100 100">
                <circle
                  cx="50"
                  cy="50"
                  r="40"
                  fill="none"
                  stroke="var(--glass-border)"
                  stroke-width="6"
                />
                <circle
                  cx="50"
                  cy="50"
                  r="40"
                  fill="none"
                  stroke="var(--accent)"
                  stroke-width="6"
                  stroke-dasharray={`${($responseTime / 2) * 251} 251`}
                  stroke-dashoffset="0"
                  transform="rotate(-90 50 50)"
                />
              </svg>
              <div class="metric-text">
                <span class="metric-number">{$responseTime.toFixed(0)}</span>
                <span class="metric-unit">s</span>
              </div>
            </div>
            <p class="metric-label">Response Time</p>
          </div>

          <div class="metric-display" in:fade={{ delay: 100 }}>
            <div class="metric-circle">
              <svg width="100" height="100" viewBox="0 0 100 100">
                <circle
                  cx="50"
                  cy="50"
                  r="40"
                  fill="none"
                  stroke="var(--glass-border)"
                  stroke-width="6"
                />
                <circle
                  cx="50"
                  cy="50"
                  r="40"
                  fill="none"
                  stroke="var(--accent)"
                  stroke-width="6"
                  stroke-dasharray={`${($concurrentUsers / 100) * 251} 251`}
                  stroke-dashoffset="0"
                  transform="rotate(-90 50 50)"
                />
              </svg>
              <div class="metric-text">
                <span class="metric-number">{$concurrentUsers.toFixed(0)}</span>
                <span class="metric-unit">+</span>
              </div>
            </div>
            <p class="metric-label">Concurrent Users</p>
          </div>

          <div class="metric-display" in:fade={{ delay: 200 }}>
            <div class="metric-circle">
              <svg width="100" height="100" viewBox="0 0 100 100">
                <circle
                  cx="50"
                  cy="50"
                  r="40"
                  fill="none"
                  stroke="var(--glass-border)"
                  stroke-width="6"
                />
                <circle
                  cx="50"
                  cy="50"
                  r="40"
                  fill="none"
                  stroke="var(--accent)"
                  stroke-width="6"
                  stroke-dasharray={`${($uptime / 99.9) * 251} 251`}
                  stroke-dashoffset="0"
                  transform="rotate(-90 50 50)"
                />
              </svg>
              <div class="metric-text">
                <span class="metric-number">{$uptime.toFixed(1)}</span>
                <span class="metric-unit">%</span>
              </div>
            </div>
            <p class="metric-label">Uptime</p>
          </div>
        </div>

        <div class="specs-grid">
          {#each performanceSpecs as spec}
            <div class="spec-card">
              <div class="spec-header">
                <h5>{spec.metric}</h5>
                <span class="spec-value">{spec.value}</span>
              </div>
              <p>{spec.description}</p>
            </div>
          {/each}
        </div>
      </div>
    {/if}
  </div>

  <div class="infrastructure-note" in:fade={{ delay: 400 }}>
    <DynamicIcon name="InfoIcon" size={20} />
    <p>Powered by cloud infrastructure with auto-scaling capabilities</p>
  </div>
</div>

<style>
  .technical-requirements {
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
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

  .tabs {
    display: flex;
    gap: 0.5rem;
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 8px;
    padding: 0.25rem;
  }

  .tab {
    flex: 1;
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.5rem;
    padding: 0.75rem 1rem;
    background: transparent;
    border: none;
    border-radius: 6px;
    color: var(--text-secondary);
    font-size: 0.9rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.3s ease;
  }

  .tab:hover {
    color: var(--text-primary);
  }

  .tab.active {
    background: var(--surface);
    color: var(--accent);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  }

  .tab-content {
    min-height: 400px;
  }

  .compatibility-section {
    display: flex;
    flex-direction: column;
    gap: 2rem;
  }

  .compatibility-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1.5rem;
  }

  .compat-category {
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 12px;
    padding: 1.25rem;
  }

  .compat-category h4 {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin: 0 0 1rem;
    font-size: 1rem;
    color: var(--text-primary);
  }

  .compat-items {
    display: flex;
    flex-direction: column;
    gap: 0.75rem;
  }

  .compat-item {
    display: flex;
    align-items: flex-start;
    gap: 0.75rem;
    padding: 0.75rem;
    background: var(--surface);
    border-radius: 8px;
    transition: all 0.3s ease;
  }

  .compat-item:hover {
    transform: translateX(4px);
  }

  .compat-item.supported {
    border-left: 3px solid var(--success);
  }

  .compat-item.limited {
    border-left: 3px solid var(--warning);
  }

  .compat-item.unsupported {
    border-left: 3px solid var(--error);
  }

  .compat-item .label {
    display: block;
    font-size: 0.9rem;
    color: var(--text-primary);
    font-weight: 500;
  }

  .compat-item .value {
    display: block;
    font-size: 0.85rem;
    color: var(--text-secondary);
    margin-top: 0.25rem;
  }

  .element-showcase {
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 12px;
    padding: 1.5rem;
    text-align: center;
  }

  .element-showcase h5 {
    margin: 0 0 1rem;
    font-size: 1.1rem;
    color: var(--text-primary);
  }

  .element-grid {
    display: flex;
    flex-wrap: wrap;
    gap: 0.75rem;
    justify-content: center;
  }

  .element-badge {
    width: 40px;
    height: 40px;
    display: flex;
    align-items: center;
    justify-content: center;
    background: var(--accent);
    color: var(--surface);
    border-radius: 8px;
    font-weight: 600;
    font-size: 1rem;
    transition: all 0.3s ease;
  }

  .element-badge:hover {
    transform: scale(1.1);
    box-shadow: 0 4px 12px rgba(var(--accent-rgb), 0.3);
  }

  .performance-section {
    display: flex;
    flex-direction: column;
    gap: 2rem;
  }

  .performance-metrics {
    display: flex;
    justify-content: space-around;
    gap: 2rem;
    flex-wrap: wrap;
  }

  .metric-display {
    text-align: center;
  }

  .metric-circle {
    position: relative;
    width: 100px;
    height: 100px;
    margin: 0 auto 0.5rem;
  }

  .metric-text {
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    display: flex;
    align-items: baseline;
    gap: 0.25rem;
  }

  .metric-number {
    font-size: 1.5rem;
    font-weight: 600;
    color: var(--text-primary);
  }

  .metric-unit {
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .metric-label {
    margin: 0;
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .specs-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 1rem;
  }

  .spec-card {
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 8px;
    padding: 1rem;
    transition: all 0.3s ease;
  }

  .spec-card:hover {
    transform: translateY(-2px);
    border-color: var(--accent-dim);
  }

  .spec-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 0.5rem;
  }

  .spec-header h5 {
    margin: 0;
    font-size: 0.9rem;
    color: var(--text-primary);
  }

  .spec-value {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--accent);
  }

  .spec-card p {
    margin: 0;
    font-size: 0.8rem;
    color: var(--text-secondary);
  }

  .infrastructure-note {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 1rem;
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 8px;
    color: var(--text-secondary);
    font-size: 0.9rem;
  }

  .infrastructure-note p {
    margin: 0;
  }

  @media (max-width: 768px) {
    .compatibility-grid {
      grid-template-columns: 1fr;
    }

    .performance-metrics {
      justify-content: center;
    }
  }
</style>
