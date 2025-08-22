<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
  import MoleculeViewer from './components/MoleculeViewer.svelte';
  import { fade, fly } from 'svelte/transition';
  import { onMount } from 'svelte';

  let activeStep = 0;
  let animationTimer: NodeJS.Timeout;

  const fragmentationSteps = [
    {
      title: 'Initial Ionization',
      description: '[M+H]+ ion formation',
      smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
      highlights: [],
    },
    {
      title: 'Bond Cleavage',
      description: 'Weakest C-N bond breaks',
      smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
      highlights: [1, 2],
    },
    {
      title: 'Fragment Formation',
      description: 'Base peak at m/z 109',
      smiles: 'CN1C=NC2=C1C(=O)N',
      highlights: [],
    },
  ];

  const fragmentationFactors = [
    {
      icon: 'BondIcon',
      title: 'Bond Strength',
      description: 'C-C > C-N > C-O',
      value: '350-250 kJ/mol',
    },
    {
      icon: 'AtomIcon',
      title: 'Charge Stability',
      description: 'Resonance stabilization',
      value: 'Aromatic > Aliphatic',
    },
    {
      icon: 'MoleculeIcon',
      title: 'Rearrangements',
      description: 'McLafferty, H-transfer',
      value: '6-member transition',
    },
  ];

  onMount(() => {
    // Auto-advance animation
    animationTimer = setInterval(() => {
      activeStep = (activeStep + 1) % fragmentationSteps.length;
    }, 3000);

    return () => {
      if (animationTimer) clearInterval(animationTimer);
    };
  });

  function selectStep(index: number) {
    activeStep = index;
    if (animationTimer) {
      clearInterval(animationTimer);
      animationTimer = setInterval(() => {
        activeStep = (activeStep + 1) % fragmentationSteps.length;
      }, 3000);
    }
  }
</script>

<div class="fragmentation-science">
  <div class="header-section">
    <DynamicIcon name="FragmentationIcon" size={40} color="var(--accent)" />
    <h3>Fragmentation Science</h3>
    <p class="subtitle">Understanding how molecules break apart in mass spectrometry</p>
  </div>

  <div class="visualization-section">
    <div class="molecule-animation">
      {#key activeStep}
        <div class="step-content" in:fly={{ y: 20, duration: 300 }}>
          <MoleculeViewer
            smiles={fragmentationSteps[activeStep].smiles}
            highlights={fragmentationSteps[activeStep].highlights}
            width={280}
            height={200}
            interactive={true}
          />
          <div class="step-info">
            <h4>{fragmentationSteps[activeStep].title}</h4>
            <p>{fragmentationSteps[activeStep].description}</p>
          </div>
        </div>
      {/key}
    </div>

    <div class="step-indicators">
      {#each fragmentationSteps as step, i}
        <button
          class="step-dot"
          class:active={activeStep === i}
          on:click={() => selectStep(i)}
          aria-label={`Step ${i + 1}: ${step.title}`}
        ></button>
      {/each}
    </div>
  </div>

  <div class="factors-grid">
    {#each fragmentationFactors as factor, i}
      <div class="factor-card hover-lift" in:fade={{ delay: i * 100, duration: 300 }}>
        <div class="factor-header">
          <DynamicIcon name={factor.icon as any} size={24} color="var(--accent)" />
          <h5>{factor.title}</h5>
        </div>
        <p class="factor-description">{factor.description}</p>
        <div class="factor-value">{factor.value}</div>
      </div>
    {/each}
  </div>

  <div class="example-showcase">
    <div class="example-header">
      <h4>Real Example: Caffeine Fragmentation</h4>
      <DynamicIcon name="ArrowRightIcon" size={20} />
    </div>
    <div class="example-details">
      <div class="mass-transition">
        <span class="mass-value">194</span>
        <span class="arrow">→</span>
        <span class="mass-value">109</span>
      </div>
      <p class="loss-info">Loss of C₄H₅NO (-85 Da)</p>
    </div>
  </div>
</div>

<style>
  /* Import only needed animations */
  .hover-lift {
    transition:
      transform var(--transition-smooth),
      box-shadow var(--transition-smooth);
  }

  .hover-lift:hover {
    transform: translateY(-2px);
    box-shadow: 0 8px 25px rgba(0, 0, 0, 0.08);
  }
  .fragmentation-science {
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

  .visualization-section {
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 12px;
    padding: 1.5rem;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 1rem;
  }

  .molecule-animation {
    position: relative;
    min-height: 280px;
    width: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .step-content {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 1rem;
  }

  .step-info {
    text-align: center;
  }

  .step-info h4 {
    margin: 0 0 0.25rem;
    font-size: 1.1rem;
    color: var(--accent);
  }

  .step-info p {
    margin: 0;
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .step-indicators {
    display: flex;
    gap: 0.5rem;
  }

  .step-dot {
    width: 8px;
    height: 8px;
    border-radius: 50%;
    border: none;
    background: var(--glass-border);
    cursor: pointer;
    transition: all 0.3s ease;
    padding: 0;
  }

  .step-dot:hover {
    background: var(--accent-dim);
    transform: scale(1.2);
  }

  .step-dot.active {
    background: var(--accent);
    width: 24px;
    border-radius: 4px;
  }

  .factors-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
  }

  .factor-card {
    background: var(--glass-bg);
    border: 1px solid var(--glass-border);
    border-radius: 8px;
    padding: 1rem;
  }

  .factor-card:hover {
    border-color: var(--accent-dim);
  }

  .factor-header {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin-bottom: 0.5rem;
  }

  .factor-header h5 {
    margin: 0;
    font-size: 1rem;
    color: var(--text-primary);
  }

  .factor-description {
    margin: 0 0 0.5rem;
    font-size: 0.85rem;
    color: var(--text-secondary);
  }

  .factor-value {
    font-size: 0.9rem;
    color: var(--accent);
    font-weight: 500;
  }

  .example-showcase {
    background: linear-gradient(135deg, var(--glass-bg), transparent);
    border: 1px solid var(--glass-border);
    border-radius: 8px;
    padding: 1.25rem;
  }

  .example-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 1rem;
  }

  .example-header h4 {
    margin: 0;
    font-size: 1.1rem;
    color: var(--text-primary);
  }

  .example-details {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
  }

  .mass-transition {
    display: flex;
    align-items: center;
    gap: 1rem;
    font-size: 1.5rem;
    font-weight: 600;
  }

  .mass-value {
    color: var(--accent);
    font-family: var(--font-mono);
  }

  .arrow {
    color: var(--text-secondary);
  }

  .loss-info {
    margin: 0;
    font-size: 0.9rem;
    color: var(--text-secondary);
    font-family: var(--font-mono);
  }

  @media (max-width: 768px) {
    .factors-grid {
      grid-template-columns: 1fr;
    }
  }
</style>
