<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;

  const supervisor = {
    name: 'Dr. Mohammed Elbadawi',
    title: 'Research Supervisor',
    department: 'Computer Science',
    institution: 'Queen Mary University of London',
  };

  const contributors = [
    { name: 'RDKit Development Team', icon: 'RDKitIcon' },
    { name: 'Scikit-learn Contributors', icon: 'ScikitLearnIcon' },
    { name: 'SvelteKit Community', icon: 'SvelteKitIcon' },
    { name: 'QMUL AI in Biosciences', icon: 'GraduationCapIcon' },
    { name: 'Mass Spectrometry Research', icon: 'MicroscopeIcon' },
  ];
</script>

<div class="acknowledgments-panel" class:compact>
  <div class="panel-content">
    <div class="hero-section">
      <div class="hero-icon">
        <DynamicIcon name="HeartIcon" size={compact ? 32 : 48} />
      </div>
      <h2 class="section-title text-gradient">Acknowledgments</h2>
      <p class="section-subtitle">
        Grateful for the guidance, tools, and community that made this project possible
      </p>
    </div>

    <div class="supervisor-card">
      <div class="supervisor-icon">
        <DynamicIcon name="GraduationCapIcon" size={compact ? 28 : 36} />
      </div>

      <div class="supervisor-info">
        <h3 class="supervisor-name">{supervisor.name}</h3>
        <div class="supervisor-title">{supervisor.title}</div>
        <div class="supervisor-meta">
          <span>{supervisor.department}</span>
          <span class="separator">•</span>
          <span>{supervisor.institution}</span>
        </div>
      </div>
    </div>

    <div class="contributors-grid">
      <h3 class="contributors-title">Open Source & Research Communities</h3>
      <div class="contributors-list">
        {#each contributors as contributor}
          <div class="contributor-item">
            <div class="contributor-icon">
              <DynamicIcon name={contributor.icon as any} size={compact ? 20 : 24} />
            </div>
            <span class="contributor-name">{contributor.name}</span>
          </div>
        {/each}
      </div>
    </div>
  </div>
</div>

<style>
  .acknowledgments-panel {
    display: flex;
    justify-content: center;
    padding: 3rem 2rem;
  }

  .acknowledgments-panel.compact {
    padding: 2rem 1rem;
  }

  .panel-content {
    display: flex;
    flex-direction: column;
    gap: 2rem;
    max-width: 900px;
    width: 100%;
  }

  .hero-section {
    text-align: center;
    margin-bottom: 1rem;
  }

  .hero-icon {
    color: var(--accent);
    margin-bottom: 1rem;
    display: flex;
    justify-content: center;
  }

  .section-title {
    font-size: 2rem;
    font-weight: 700;
    margin: 0 0 0.5rem 0;
    letter-spacing: -0.02em;
  }

  .section-subtitle {
    font-size: 1rem;
    color: var(--text-secondary);
    margin: 0;
    max-width: 600px;
    margin-left: auto;
    margin-right: auto;
    line-height: 1.5;
  }

  .text-gradient {
    background: linear-gradient(to right, var(--accent), var(--accent-secondary));
    -webkit-background-clip: text;
    background-clip: text;
    -webkit-text-fill-color: transparent;
    color: transparent;
  }

  .supervisor-card {
    display: flex;
    gap: 1.5rem;
    align-items: center;
    background: var(--surface-glass);
    border-radius: var(--radius-lg);
    padding: 2rem;
    border: 1px solid rgba(255, 255, 255, 0.2);
    transition:
      transform 0.3s ease,
      box-shadow 0.3s ease;
    position: relative;
    overflow: hidden;
  }

  .supervisor-card:hover {
    transform: translateY(-5px);
    box-shadow: 0 15px 30px rgba(0, 0, 0, 0.1);
  }

  .supervisor-icon {
    width: 60px;
    height: 60px;
    background: var(--accent-soft);
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    color: var(--accent);
    flex-shrink: 0;
  }

  .supervisor-info {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
  }

  .supervisor-name {
    font-size: 1.4rem;
    font-weight: 600;
    margin: 0;
    color: var(--text-primary);
  }

  .supervisor-title {
    font-size: 1rem;
    color: var(--accent);
    font-weight: 500;
  }

  .supervisor-meta {
    font-size: 0.9rem;
    color: var(--text-secondary);
    display: flex;
    align-items: center;
    gap: 0.5rem;
  }

  .separator {
    opacity: 0.6;
  }

  .contributors-grid {
    background: var(--surface-glass);
    border-radius: var(--radius-lg);
    padding: 2rem;
    border: 1px solid rgba(255, 255, 255, 0.2);
  }

  .contributors-title {
    font-size: 1.1rem;
    font-weight: 600;
    margin: 0 0 1.5rem 0;
    color: var(--text-primary);
    text-align: center;
  }

  .contributors-list {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
  }

  .contributor-item {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 0.75rem 1rem;
    background: rgba(255, 255, 255, 0.05);
    border-radius: var(--radius-md);
    transition: all 0.3s ease;
    border: 1px solid transparent;
  }

  .contributor-item:hover {
    background: rgba(120, 121, 255, 0.1);
    border-color: rgba(120, 121, 255, 0.3);
    transform: translateY(-2px);
  }

  .contributor-icon {
    color: var(--accent);
    display: flex;
    align-items: center;
    justify-content: center;
    flex-shrink: 0;
  }

  .contributor-name {
    font-size: 0.9rem;
    color: var(--text-secondary);
    font-weight: 500;
  }

  .contributor-item:hover .contributor-name {
    color: var(--accent);
  }

  @media (max-width: 768px) {
    .acknowledgments-panel {
      padding: 2rem 1rem;
    }

    .section-title {
      font-size: 1.6rem;
    }

    .section-subtitle {
      font-size: 0.9rem;
    }

    .supervisor-card {
      flex-direction: column;
      text-align: center;
      gap: 1rem;
      padding: 1.5rem;
    }

    .supervisor-name {
      font-size: 1.2rem;
    }

    .supervisor-meta {
      flex-direction: column;
      gap: 0.25rem;
    }

    .separator {
      display: none;
    }

    .contributors-list {
      grid-template-columns: 1fr;
    }

    .contributor-item {
      padding: 1rem;
    }
  }

  @media (max-width: 480px) {
    .contributors-grid {
      padding: 1.5rem;
    }

    .contributors-title {
      font-size: 1rem;
    }
  }
</style>
