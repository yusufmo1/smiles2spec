<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;

  const socialLinks = [
    { icon: 'GithubIcon', url: 'https://github.com/yusufmo1', label: 'GitHub' },
    {
      icon: 'LinkedinIcon',
      url: 'https://www.linkedin.com/in/yusuf-mohammed1/',
      label: 'LinkedIn',
    },
    { icon: 'MailIcon', url: 'mailto:contact@email.com', label: 'Email' },
  ];
</script>

<div class="hero-panel" class:compact>
  <div class="hero-content">
    <div class="photo-section">
      <div class="photo-frame">
        <img src="/assets/images/yusuf-photo.png" alt="Yusuf Mohammed" />
        <div class="photo-glow"></div>
      </div>
      <div class="floating-particles">
        {#each Array(6) as _, i}
          <div class="particle" style="--delay: {i * 0.5}s"></div>
        {/each}
      </div>
    </div>

    <div class="info-section">
      <h1 class="developer-name text-gradient">Yusuf Mohammed</h1>
      <h2 class="developer-title">MSc AI in the Biosciences</h2>
      <p class="developer-bio">
        Bridging biochemistry and machine learning to democratize mass spectrometry analysis.
        Passionate about making complex scientific tools accessible to researchers worldwide.
      </p>

      <div class="social-links">
        {#each socialLinks as link}
          <a
            href={link.url}
            target="_blank"
            rel="noopener noreferrer"
            class="social-link"
            aria-label={link.label}
          >
            <DynamicIcon name={link.icon as any} size={24} />
          </a>
        {/each}
      </div>
    </div>
  </div>
</div>

<style>
  .hero-panel {
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 320px;
    padding: 1.5rem;
  }

  .hero-content {
    display: grid;
    grid-template-columns: 40% 60%;
    gap: 1.5rem;
    align-items: center;
    max-width: 1000px;
    width: 100%;
  }

  .photo-section {
    position: relative;
  }

  .photo-frame {
    position: relative;
    width: 200px;
    height: 200px;
    margin: 0 auto;
    border-radius: 50%;
    overflow: hidden;
    background: var(--surface-glass);
    border: 3px solid rgba(120, 121, 255, 0.3);
    box-shadow: 0 20px 40px rgba(120, 121, 255, 0.2);
  }

  .photo-frame img {
    width: 100%;
    height: 100%;
    object-fit: cover;
  }

  .photo-glow {
    position: absolute;
    inset: -50%;
    background: radial-gradient(circle, var(--accent-soft) 0%, transparent 70%);
    animation: rotate 20s linear infinite;
    z-index: -1;
  }

  .floating-particles {
    position: absolute;
    inset: 0;
    pointer-events: none;
  }

  .particle {
    position: absolute;
    width: 4px;
    height: 4px;
    background: var(--accent);
    border-radius: 50%;
    opacity: 0.6;
    animation: float 6s ease-in-out infinite;
    animation-delay: var(--delay);
  }

  .particle:nth-child(1) {
    top: 10%;
    left: 20%;
  }
  .particle:nth-child(2) {
    top: 20%;
    right: 10%;
  }
  .particle:nth-child(3) {
    bottom: 30%;
    left: 10%;
  }
  .particle:nth-child(4) {
    bottom: 10%;
    right: 20%;
  }
  .particle:nth-child(5) {
    top: 50%;
    left: 5%;
  }
  .particle:nth-child(6) {
    top: 60%;
    right: 5%;
  }

  .developer-name {
    font-size: 2.5rem;
    font-weight: 700;
    margin: 0 0 0.25rem 0;
    letter-spacing: -0.02em;
  }

  .developer-title {
    font-size: 1.25rem;
    color: var(--accent);
    margin: 0 0 0.75rem 0;
    font-weight: 500;
  }

  .developer-bio {
    font-size: 1rem;
    line-height: 1.5;
    color: var(--text-secondary);
    margin: 0 0 1rem 0;
  }

  .social-links {
    display: flex;
    gap: 0.75rem;
  }

  .social-link {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 40px;
    height: 40px;
    background: var(--surface-glass);
    border: 1px solid rgba(120, 121, 255, 0.2);
    border-radius: 50%;
    color: var(--accent);
    transition: all 0.3s ease;
  }

  .social-link:hover {
    background: var(--accent);
    color: white;
    transform: translateY(-4px);
    box-shadow: 0 10px 20px rgba(120, 121, 255, 0.3);
  }

  .text-gradient {
    background: linear-gradient(to right, var(--accent), var(--accent-secondary));
    -webkit-background-clip: text;
    background-clip: text;
    -webkit-text-fill-color: transparent;
    color: transparent;
  }

  @keyframes rotate {
    to {
      transform: rotate(360deg);
    }
  }

  @keyframes float {
    0%,
    100% {
      transform: translateY(0) scale(1);
    }
    50% {
      transform: translateY(-20px) scale(1.2);
    }
  }

  @media (max-width: 768px) {
    .hero-panel {
      padding: 1rem;
    }

    .hero-content {
      grid-template-columns: 1fr;
      text-align: center;
      gap: 1rem;
    }

    .photo-frame {
      width: 180px;
      height: 180px;
    }

    .developer-name {
      font-size: 2rem;
    }

    .developer-title {
      font-size: 1.1rem;
      margin-bottom: 0.5rem;
    }

    .developer-bio {
      font-size: 0.9rem;
      margin-bottom: 0.75rem;
    }

    .social-links {
      justify-content: center;
    }
  }
</style>
