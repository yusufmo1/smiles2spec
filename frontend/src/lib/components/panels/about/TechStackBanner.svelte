<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;

  const frontendTech = [
    { icon: 'NPMIcon', name: 'NPM' },
    { icon: 'SvelteKitIcon', name: 'SvelteKit' },
    { icon: 'TypeScriptIcon', name: 'TypeScript' },
    { icon: 'CSSIcon', name: 'CSS3' },
    { icon: 'PlotlyJSIcon', name: 'Plotly.js' },
  ];

  const backendTech = [
    { icon: 'PythonIcon', name: 'Python' },
    { icon: 'FlaskIcon', name: 'Flask' },
    { icon: 'RDKitIcon', name: 'RDKit' },
    { icon: 'ScikitLearnIcon', name: 'Scikit-learn' },
    { icon: 'MatplotlibIcon', name: 'Matplotlib' },
    { icon: 'PyTorchIcon', name: 'PyTorch' },
    { icon: 'TensorFlowIcon', name: 'TensorFlow' },
    { icon: 'DockerIcon', name: 'Docker' },
    { icon: 'GithubIcon', name: 'GitHub Actions' },
  ];
</script>

<div class="tech-banners-container" class:compact>
  <!-- Frontend Panel -->
  <div class="tech-banner frontend-banner">
    <h3 class="panel-title">Frontend Stack</h3>
    <div class="scroll-container" role="region" aria-label="frontend tech stack carousel">
      <div class="tech-track scroll-right">
        {#each [0, 1, 2] as _}
          <div class="category-section">
            {#each frontendTech as tech}
              <div class="tech-item">
                <DynamicIcon name={tech.icon as any} size={40} />
                <span class="tech-name">{tech.name}</span>
              </div>
            {/each}
          </div>
        {/each}
      </div>
    </div>
  </div>

  <!-- Backend Panel -->
  <div class="tech-banner backend-banner">
    <h3 class="panel-title">Backend Stack</h3>
    <div class="scroll-container" role="region" aria-label="backend tech stack carousel">
      <div class="tech-track scroll-left">
        {#each [0, 1, 2] as _}
          <div class="category-section">
            {#each backendTech as tech}
              <div class="tech-item">
                <DynamicIcon name={tech.icon as any} size={40} />
                <span class="tech-name">{tech.name}</span>
              </div>
            {/each}
          </div>
        {/each}
      </div>
    </div>
  </div>
</div>

<style>
  .tech-banners-container {
    display: flex;
    flex-direction: column;
    gap: 0.75rem;
  }

  .tech-banner {
    position: relative;
    background: var(--surface-glass);
    border-radius: var(--radius-pill);
    padding: 0.75rem 0;
    overflow: hidden;
  }

  .panel-title {
    text-align: center;
    margin-bottom: 0.25rem;
    font-size: 1rem;
    font-weight: 600;
    color: var(--text-primary);
  }

  .scroll-container {
    overflow-x: hidden;
    scroll-behavior: auto;
  }

  .tech-track {
    display: flex;
    gap: 2rem;
    padding: 0 1rem;
    width: fit-content;
    animation-duration: 30s;
    animation-timing-function: linear;
    animation-iteration-count: infinite;
  }

  .scroll-right {
    animation-name: scroll-right;
  }

  .scroll-left {
    animation-name: scroll-left;
  }

  @keyframes scroll-right {
    0% {
      transform: translateX(0);
    }
    100% {
      transform: translateX(-33.333%);
    }
  }

  @keyframes scroll-left {
    0% {
      transform: translateX(-33.333%);
    }
    100% {
      transform: translateX(0);
    }
  }

  .category-section {
    display: flex;
    align-items: center;
    gap: 1.5rem;
  }

  .tech-item {
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0.3rem;
    padding: 0.5rem;
  }

  .tech-name {
    font-size: 0.8rem;
    color: var(--text-secondary);
    opacity: 0.7;
  }

  /* Add a glowing effect to icons */
  .tech-item :global(svg) {
    filter: drop-shadow(0 0 5px rgba(120, 121, 255, 0.3));
  }

  /* Frontend banner glow effect */
  .frontend-banner .tech-item :global(svg) {
    filter: drop-shadow(0 0 5px rgba(120, 121, 255, 0.4));
  }

  /* Backend banner glow effect */
  .backend-banner .tech-item :global(svg) {
    filter: drop-shadow(0 0 5px rgba(72, 187, 120, 0.4));
  }

  @media (max-width: 768px) {
    .tech-banner {
      padding: 0.5rem 0;
    }

    .category-section {
      gap: 1rem;
    }

    .tech-item {
      padding: 0.4rem;
    }

    .panel-title {
      font-size: 0.9rem;
      margin-bottom: 0.2rem;
    }
  }
</style>
