<script lang="ts">
  import { DynamicIcon, type IconName } from '$lib/components/icons';

  interface QuickLink {
    href: string;
    icon: IconName;
    title: string;
    description: string;
    color: string;
  }

  const quickLinks: QuickLink[] = [
    {
      href: '/spectral-simulation',
      icon: 'SpectrumIcon',
      title: 'Spectral Simulation',
      description:
        'Interactive prediction interface with real-time visualization, bulk processing, and multiple export formats',
      color: 'var(--accent)',
    },
    {
      href: '/chat-with-spectrum',
      icon: 'ChatWithSpectrumIcon',
      title: 'Chat with Spectrum',
      description:
        'AI-powered spectral analysis assistant with context-aware interpretation and educational guidance',
      color: '#34c759',
    },
    {
      href: '/how-it-works',
      icon: 'BookIcon',
      title: 'How It Works',
      description:
        'Deep dive into machine learning methodology, feature engineering, and validation approaches',
      color: 'var(--accent-secondary)',
    },
    {
      href: '/about',
      icon: 'UserIcon',
      title: 'About & Source',
      description:
        'Project development story, open-source contributions, and technical documentation access',
      color: '#ff6b35',
    },
  ];
</script>

<section class="quick-links-section">
  <div class="quick-links-header">
    <div class="title-pill">
      <div class="pill-title">Explore the Platform</div>
      <div class="pill-subtitle">
        Comprehensive tools and resources for mass spectrum prediction, analysis, and research
      </div>
    </div>
  </div>

  <div class="quick-links-grid">
    {#each quickLinks as link}
      <a href={link.href} class="quick-link-card glass-card" data-sveltekit-preload-data="hover">
        <div class="quick-link-icon" style="color: {link.color}">
          <DynamicIcon name={link.icon} size={28} color={link.color} />
        </div>
        <h3 class="quick-link-title">{link.title}</h3>
        <p class="quick-link-description">{link.description}</p>
        <div class="quick-link-arrow">
          <DynamicIcon name="ArrowRightIcon" size={16} color="var(--text-secondary)" />
        </div>
      </a>
    {/each}
  </div>

  <div class="footer-info">
    <p class="footer-text">
      Open source • Built with <DynamicIcon name="HeartIcon" size={14} color="#ff453a" />
      for the scientific community
    </p>
  </div>
</section>

<style>
  .quick-links-section {
    padding: 2rem 1rem 1rem;
    max-width: 1000px;
    margin: 0 auto;
  }

  .quick-links-header {
    text-align: center;
    margin-bottom: 3rem;
  }

  .title-pill {
    display: inline-flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
    padding: 2rem 3rem;
    background: linear-gradient(
      135deg,
      rgba(120, 121, 255, 0.08) 0%,
      rgba(164, 93, 230, 0.08) 100%
    );
    border: 1px solid rgba(120, 121, 255, 0.2);
    border-radius: 32px;
    margin: 0 0 2rem 0;
    backdrop-filter: blur(20px);
    -webkit-backdrop-filter: blur(20px);
    transition: all var(--transition-smooth);
    position: relative;
    overflow: hidden;
    max-width: 700px;
    box-shadow:
      0 8px 32px rgba(120, 121, 255, 0.1),
      inset 0 1px 0 rgba(255, 255, 255, 0.1);
  }

  .title-pill::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: linear-gradient(
      135deg,
      rgba(120, 121, 255, 0.05) 0%,
      rgba(164, 93, 230, 0.05) 100%
    );
    border-radius: inherit;
    opacity: 0;
    transition: opacity var(--transition-smooth);
  }

  .title-pill:hover {
    transform: translateY(-2px);
    box-shadow:
      0 12px 40px rgba(120, 121, 255, 0.15),
      inset 0 1px 0 rgba(255, 255, 255, 0.15);
  }

  .title-pill:hover::before {
    opacity: 1;
  }

  .pill-title {
    font-size: clamp(1.5rem, 4vw, 2rem);
    font-weight: 700;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    -webkit-background-clip: text;
    background-clip: text;
    -webkit-text-fill-color: transparent;
    color: transparent;
    margin-bottom: 0.75rem;
    line-height: 1.2;
    position: relative;
    z-index: 1;
  }

  .pill-subtitle {
    font-size: clamp(0.95rem, 2.2vw, 1.1rem);
    color: var(--text-secondary);
    line-height: 1.5;
    font-weight: 400;
    position: relative;
    z-index: 1;
  }

  .quick-links-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
    gap: 1.5rem;
    margin-bottom: 3rem;
  }

  .quick-link-card {
    display: block;
    padding: 2rem;
    border-radius: 20px;
    background: rgba(255, 255, 255, 0.06);
    backdrop-filter: blur(15px);
    -webkit-backdrop-filter: blur(15px);
    border: 1px solid rgba(255, 255, 255, 0.08);
    text-decoration: none;
    transition: all var(--transition-smooth);
    position: relative;
    overflow: hidden;
    cursor: pointer;
  }

  .quick-link-card::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: linear-gradient(
      135deg,
      rgba(255, 255, 255, 0.08) 0%,
      rgba(255, 255, 255, 0.04) 100%
    );
    border-radius: inherit;
    opacity: 0;
    transition: opacity var(--transition-smooth);
  }

  .quick-link-card:hover {
    transform: translateY(-2px);
    box-shadow: 0 16px 32px rgba(0, 0, 0, 0.08);
    border-color: rgba(255, 255, 255, 0.15);
  }

  .quick-link-card:hover::before {
    opacity: 1;
  }

  .quick-link-card:hover .quick-link-arrow {
    transform: translateX(4px);
  }

  .quick-link-icon {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 56px;
    height: 56px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 14px;
    margin-bottom: 1.5rem;
    position: relative;
    z-index: 1;
  }

  .quick-link-title {
    font-size: 1.125rem;
    font-weight: 600;
    color: var(--text-primary);
    margin: 0 0 0.75rem 0;
    position: relative;
    z-index: 1;
  }

  .quick-link-description {
    font-size: 0.9rem;
    color: var(--text-secondary);
    line-height: 1.5;
    margin: 0 0 1rem 0;
    position: relative;
    z-index: 1;
  }

  .quick-link-arrow {
    display: flex;
    align-items: center;
    justify-content: flex-end;
    transition: transform var(--transition-smooth);
    position: relative;
    z-index: 1;
  }

  .footer-info {
    text-align: center;
    padding-top: 2rem;
    border-top: 1px solid rgba(255, 255, 255, 0.1);
  }

  .footer-text {
    font-size: 0.875rem;
    color: var(--text-secondary);
    margin: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.25rem;
    flex-wrap: wrap;
  }

  /* Responsive adjustments */
  @media (max-width: 768px) {
    .quick-links-section {
      padding: 2rem 1rem 1rem;
    }

    .title-pill {
      padding: 1.5rem 2rem;
      max-width: 90%;
      margin: 0 auto 2rem auto;
    }

    .quick-links-grid {
      grid-template-columns: 1fr;
      gap: 1rem;
    }

    .quick-link-card {
      padding: 1.5rem;
    }
  }

  @media (max-width: 480px) {
    .title-pill {
      padding: 1.25rem 1.5rem;
      border-radius: 24px;
    }

    .pill-title {
      margin-bottom: 0.5rem;
    }

    .quick-link-card {
      padding: 1.25rem;
    }

    .quick-link-icon {
      width: 48px;
      height: 48px;
      margin-bottom: 1rem;
    }

    .footer-text {
      font-size: 0.8rem;
    }
  }
</style>
