<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;

  const contactMethods = [
    {
      id: 'github',
      title: 'GitHub',
      description: 'Check out the source code, contribute, or raise issues',
      icon: 'GithubIcon',
      link: 'https://github.com/yusufmo1',
      themeClass: 'github-theme',
    },
    {
      id: 'linkedin',
      title: 'LinkedIn',
      description: 'Connect professionally for discussions or collaborations',
      icon: 'LinkedinIcon',
      link: 'https://www.linkedin.com/in/yusuf-mohammed1/',
      themeClass: 'linkedin-theme',
    },
    {
      id: 'email',
      title: 'Email',
      description: 'Direct contact for inquiries or research discussions',
      icon: 'MailIcon',
      link: 'mailto:contact@email.com',
      themeClass: 'email-theme',
    },
  ];
</script>

<div class="contact-connect-panel" class:compact>
  <div class="contact-cards">
    {#each contactMethods as method}
      <a
        href={method.link}
        target="_blank"
        rel="noopener noreferrer"
        class="contact-card {method.themeClass}"
      >
        <div class="card-content">
          <div class="card-icon">
            <DynamicIcon name={method.icon as any} size={36} />
          </div>

          <h3 class="card-title">{method.title}</h3>
          <p class="card-description">{method.description}</p>

          <div class="card-action">
            <span>Connect</span>
            <div class="arrow-icon">→</div>
          </div>
        </div>

        <div class="card-decoration"></div>
      </a>
    {/each}
  </div>
</div>

<style>
  .contact-connect-panel {
    padding: 3rem 2rem;
    display: flex;
    justify-content: center;
  }

  .contact-connect-panel.compact {
    padding: 2rem 1rem;
  }

  .contact-cards {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
    gap: 2rem;
    width: 100%;
    max-width: 1200px;
  }

  .contact-card {
    position: relative;
    min-height: 280px;
    border-radius: var(--radius-lg);
    padding: 2rem;
    display: flex;
    flex-direction: column;
    text-decoration: none;
    overflow: hidden;
    transition: all 0.3s ease;
    z-index: 1;
  }

  .card-content {
    position: relative;
    z-index: 2;
    display: flex;
    flex-direction: column;
    height: 100%;
  }

  .card-decoration {
    position: absolute;
    z-index: 1;
  }

  /* GitHub Theme */
  .github-theme {
    background: #0d1117;
    color: #e6edf3;
    border: 1px solid #30363d;
  }

  .github-theme .card-icon {
    color: #f0f6fc;
  }

  .github-theme .card-title {
    color: #f0f6fc;
  }

  .github-theme .card-description {
    color: #8b949e;
  }

  .github-theme .card-action {
    color: #58a6ff;
  }

  .github-theme .card-decoration {
    bottom: 0;
    right: 0;
    width: 180px;
    height: 180px;
    background-image: radial-gradient(circle, rgba(88, 166, 255, 0.15) 0%, transparent 70%);
  }

  .github-theme:hover {
    transform: translateY(-5px);
    box-shadow: 0 10px 25px rgba(0, 0, 0, 0.2);
    border-color: #58a6ff;
  }

  /* LinkedIn Theme */
  .linkedin-theme {
    background: linear-gradient(135deg, #0077b5, #0a66c2);
    color: white;
  }

  .linkedin-theme .card-icon,
  .linkedin-theme .card-title {
    color: white;
  }

  .linkedin-theme .card-description {
    color: rgba(255, 255, 255, 0.85);
  }

  .linkedin-theme .card-action {
    color: white;
    font-weight: 600;
  }

  .linkedin-theme .card-decoration {
    top: 0;
    right: 0;
    width: 150px;
    height: 150px;
    background-image: linear-gradient(45deg, transparent 0%, rgba(255, 255, 255, 0.1) 100%);
    border-radius: 0 0 0 150px;
  }

  .linkedin-theme:hover {
    transform: translateY(-5px);
    box-shadow: 0 10px 25px rgba(0, 119, 181, 0.4);
  }

  /* Email Theme */
  .email-theme {
    background: var(--surface-glass);
    color: var(--text-primary);
    border: 1px solid rgba(255, 255, 255, 0.2);
  }

  .email-theme .card-icon {
    color: var(--accent);
  }

  .email-theme .card-title {
    color: var(--text-primary);
  }

  .email-theme .card-description {
    color: var(--text-secondary);
  }

  .email-theme .card-action {
    color: var(--accent);
  }

  .email-theme .card-decoration {
    inset: 0;
    background: linear-gradient(
      135deg,
      rgba(var(--accent-rgb, 120, 121, 255), 0.05) 0%,
      rgba(var(--accent-rgb, 120, 121, 255), 0.1) 100%
    );
  }

  .email-theme:hover {
    transform: translateY(-5px);
    box-shadow: 0 10px 25px rgba(0, 0, 0, 0.1);
    border-color: var(--accent);
  }

  /* Common Card Styles */
  .card-icon {
    margin-bottom: 1.5rem;
  }

  .card-title {
    font-size: 1.4rem;
    font-weight: 600;
    margin: 0 0 1rem 0;
  }

  .card-description {
    font-size: 0.95rem;
    line-height: 1.6;
    margin: 0 0 2rem 0;
    flex-grow: 1;
  }

  .card-action {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    margin-top: auto;
  }

  .arrow-icon {
    transition: transform 0.3s ease;
  }

  .contact-card:hover .arrow-icon {
    transform: translateX(5px);
  }

  @media (max-width: 768px) {
    .contact-cards {
      grid-template-columns: 1fr;
    }

    .contact-card {
      min-height: 200px;
      padding: 1.5rem;
    }
  }
</style>
