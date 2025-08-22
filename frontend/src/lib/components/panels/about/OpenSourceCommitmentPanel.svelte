<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;

  const repoOwner = 'yusufmo1';
  const repoName = 'smiles2spec';
  const license = 'MIT License';

  const openSourceValues = [
    {
      icon: 'BookIcon',
      title: 'Transparent Science',
      description:
        'All code, models, and methodologies are openly available for review and validation',
    },
    {
      icon: 'UserIcon',
      title: 'Community Driven',
      description:
        'Built with and for the scientific community, encouraging collaboration and contribution',
    },
    {
      icon: 'HeartIcon',
      title: 'Free to Use',
      description: 'No licensing fees or restrictions - accessible to researchers worldwide',
    },
    {
      icon: 'RefreshIcon',
      title: 'Continuous Improvement',
      description: 'Open development process allows for rapid iteration and community feedback',
    },
  ];
</script>

<div class="open-source-panel" class:compact>
  <div class="panel-content">
    <div class="header-section">
      <div class="panel-header">
        <DynamicIcon name="HeartIcon" size={28} />
        <h3>Open Source Commitment</h3>
      </div>
      <p class="panel-description">
        SMILES2SPEC is built on the principles of open science and collaborative research, making
        advanced mass spectrometry prediction accessible to everyone.
      </p>
      <div class="license-badge">{license}</div>
    </div>

    <div class="values-grid">
      {#each openSourceValues as value}
        <div class="value-card">
          <div class="value-icon">
            <DynamicIcon name={value.icon as any} size={20} />
          </div>
          <div class="value-content">
            <h4 class="value-title">{value.title}</h4>
            <p class="value-description">{value.description}</p>
          </div>
        </div>
      {/each}
    </div>

    <a
      href="https://github.com/{repoOwner}/{repoName}"
      target="_blank"
      rel="noopener noreferrer"
      class="contribute-button"
    >
      <DynamicIcon name="GithubIcon" size={18} />
      <span>Explore on GitHub</span>
    </a>
  </div>
</div>

<style>
  .open-source-panel {
    padding: 3rem 2rem;
    display: flex;
    justify-content: center;
  }

  .open-source-panel.compact {
    padding: 2rem 1rem;
  }

  .panel-content {
    display: flex;
    flex-direction: column;
    gap: 2rem;
    max-width: 800px;
    width: 100%;
  }

  .header-section {
    background: var(--surface-glass);
    border-radius: var(--radius-pill);
    padding: 2rem;
    border: 1px solid rgba(255, 255, 255, 0.2);
    position: relative;
    text-align: center;
  }

  .panel-header {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.75rem;
    margin-bottom: 1rem;
    color: var(--text-primary);
  }

  .panel-header h3 {
    font-size: 1.4rem;
    font-weight: 600;
    margin: 0;
  }

  .panel-description {
    color: var(--text-secondary);
    line-height: 1.6;
    margin: 0 0 1.5rem 0;
    max-width: 600px;
    margin-left: auto;
    margin-right: auto;
  }

  .license-badge {
    background: var(--accent-soft);
    color: var(--accent);
    font-size: 0.8rem;
    font-weight: 600;
    padding: 0.4rem 1rem;
    border-radius: var(--radius-pill);
    display: inline-block;
  }

  .values-grid {
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 1rem;
  }

  .value-card {
    background: rgba(255, 255, 255, 0.6);
    border-radius: var(--radius-lg);
    padding: 1.25rem;
    border: 1px solid rgba(255, 255, 255, 0.8);
    transition: all 0.3s ease;
    display: flex;
    gap: 1rem;
    align-items: flex-start;
  }

  .value-card:hover {
    transform: translateY(-2px);
    box-shadow: 0 6px 20px rgba(0, 0, 0, 0.1);
    background: rgba(255, 255, 255, 0.8);
  }

  .value-icon {
    background: var(--accent-soft);
    color: var(--accent);
    padding: 0.75rem;
    border-radius: var(--radius-md);
    display: flex;
    align-items: center;
    justify-content: center;
    flex-shrink: 0;
  }

  .value-content {
    flex: 1;
  }

  .value-title {
    font-size: 1rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 0.5rem 0;
  }

  .value-description {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
    margin: 0;
  }

  .contribute-button {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.75rem;
    background: var(--accent);
    color: white;
    padding: 1rem 1.5rem;
    border-radius: var(--radius-pill);
    font-weight: 600;
    text-align: center;
    text-decoration: none;
    transition: all 0.3s ease;
    border: none;
    cursor: pointer;
    width: fit-content;
    align-self: center;
    box-shadow: 0 4px 10px rgba(120, 121, 255, 0.3);
  }

  .contribute-button:hover {
    background: var(--accent-hover, #6162cc);
    box-shadow: 0 6px 15px rgba(120, 121, 255, 0.4);
    transform: translateY(-2px);
  }

  @media (max-width: 768px) {
    .values-grid {
      grid-template-columns: 1fr;
      gap: 0.75rem;
    }

    .value-card {
      flex-direction: column;
      text-align: center;
      gap: 0.75rem;
    }

    .panel-header h3 {
      font-size: 1.2rem;
    }

    .header-section {
      padding: 1.5rem;
    }
  }
</style>
