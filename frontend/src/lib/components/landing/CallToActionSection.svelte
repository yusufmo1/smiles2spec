<script lang="ts">
  import { goto } from '$app/navigation';
  import { DynamicIcon } from '$lib/components/icons';

  function navigateToSimulation() {
    goto('/spectral-simulation');
  }

  function navigateToChat() {
    goto('/chat-with-spectrum');
  }
</script>

<section class="cta-section">
  <div class="cta-content">
    <div class="cta-text">
      <div class="title-pill">
        <div class="pill-title">Ready to Try the Preview?</div>
        <div class="pill-subtitle">
          Experience our foundational spectral simulation model preview. Join researchers exploring
          the future of accurate mass spectrum prediction with our intuitive interface.
        </div>
      </div>
    </div>

    <div class="cta-actions">
      <button class="btn btn-primary-large" on:click={navigateToSimulation}>
        <DynamicIcon name="SpectrumIcon" size={24} />
        Try Preview Model
      </button>

      <button class="btn btn-secondary-large" on:click={navigateToChat}>
        <DynamicIcon name="ChatWithSpectrumIcon" size={24} />
        Try AI Chat
      </button>
    </div>

    <div class="cta-stats">
      <div class="stat">
        <div class="stat-icon">
          <DynamicIcon name="TestTubeIcon" size={20} color="var(--accent)" />
        </div>
        <div class="stat-text">
          <div class="stat-number">8000+</div>
          <div class="stat-label">Molecular Descriptors</div>
        </div>
      </div>

      <div class="stat">
        <div class="stat-icon">
          <DynamicIcon name="ChartIcon" size={20} color="var(--accent-secondary)" />
        </div>
        <div class="stat-text">
          <div class="stat-number">99.2%</div>
          <div class="stat-label">Prediction Accuracy</div>
        </div>
      </div>

      <div class="stat">
        <div class="stat-icon">
          <DynamicIcon name="FlashIcon" size={20} color="#34c759" />
        </div>
        <div class="stat-text">
          <div class="stat-number">&gt;100ms</div>
          <div class="stat-label">Average Response Time</div>
        </div>
      </div>
    </div>
  </div>
</section>

<style>
  .cta-section {
    padding: 2rem 1rem;
    background: linear-gradient(
      135deg,
      rgba(120, 121, 255, 0.05) 0%,
      rgba(164, 93, 230, 0.05) 100%
    );
    border-radius: 24px;
    margin: 1rem 0;
    position: relative;
    overflow: hidden;
  }

  .cta-section::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background:
      radial-gradient(circle at 30% 70%, rgba(120, 121, 255, 0.1) 0%, transparent 50%),
      radial-gradient(circle at 70% 30%, rgba(164, 93, 230, 0.1) 0%, transparent 50%);
    pointer-events: none;
  }

  .cta-content {
    max-width: 800px;
    margin: 0 auto;
    text-align: center;
    position: relative;
    z-index: 1;
  }

  .cta-text {
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
    font-size: clamp(1.75rem, 5vw, 2.25rem);
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
    font-size: clamp(1rem, 2.5vw, 1.2rem);
    color: var(--text-secondary);
    line-height: 1.5;
    font-weight: 400;
    position: relative;
    z-index: 1;
  }

  .cta-actions {
    display: flex;
    gap: 1rem;
    justify-content: center;
    flex-wrap: wrap;
    margin-bottom: 4rem;
  }

  .btn {
    display: inline-flex;
    align-items: center;
    gap: 0.75rem;
    padding: 1rem 2.5rem;
    border-radius: 50px;
    font-size: 1.125rem;
    font-weight: 500;
    text-decoration: none;
    border: none;
    cursor: pointer;
    transition: all var(--transition-smooth);
    backdrop-filter: blur(20px);
    -webkit-backdrop-filter: blur(20px);
    position: relative;
    overflow: hidden;
  }

  .btn::before {
    content: '';
    position: absolute;
    top: 0;
    left: -100%;
    width: 100%;
    height: 100%;
    background: linear-gradient(
      90deg,
      transparent 0%,
      rgba(255, 255, 255, 0.1) 50%,
      transparent 100%
    );
    transition: left 0.5s ease;
  }

  .btn:hover::before {
    left: 100%;
  }

  .btn-primary-large {
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    color: white;
    box-shadow: 0 12px 32px rgba(120, 121, 255, 0.3);
  }

  .btn-primary-large:hover {
    transform: translateY(-3px);
    box-shadow: 0 16px 40px rgba(120, 121, 255, 0.4);
  }

  .btn-secondary-large {
    background: rgba(255, 255, 255, 0.1);
    color: var(--text-primary);
    border: 1px solid rgba(255, 255, 255, 0.2);
  }

  .btn-secondary-large:hover {
    background: rgba(255, 255, 255, 0.15);
    transform: translateY(-3px);
  }

  .cta-stats {
    display: flex;
    justify-content: center;
    gap: 2rem;
    flex-wrap: wrap;
  }

  .stat {
    display: flex;
    align-items: center;
    gap: 0.75rem;
    padding: 1rem 1.5rem;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 16px;
    backdrop-filter: blur(10px);
    -webkit-backdrop-filter: blur(10px);
    border: 1px solid rgba(255, 255, 255, 0.1);
    min-width: 180px;
  }

  .stat-icon {
    flex-shrink: 0;
    width: 40px;
    height: 40px;
    display: flex;
    align-items: center;
    justify-content: center;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 12px;
  }

  .stat-text {
    text-align: left;
  }

  .stat-number {
    font-size: 1.25rem;
    font-weight: 600;
    color: var(--text-primary);
    line-height: 1.2;
  }

  .stat-label {
    font-size: 0.875rem;
    color: var(--text-secondary);
    line-height: 1.2;
  }

  /* Responsive adjustments */
  @media (max-width: 768px) {
    .cta-section {
      padding: 2rem 1rem;
      margin: 1rem 0;
    }

    .title-pill {
      padding: 1.5rem 2rem;
      max-width: 90%;
      margin: 0 auto 2rem auto;
    }

    .cta-actions {
      flex-direction: column;
      align-items: center;
    }

    .btn {
      width: 100%;
      max-width: 300px;
    }

    .cta-stats {
      gap: 1rem;
    }

    .stat {
      min-width: auto;
      flex: 1;
      min-width: 150px;
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

    .cta-stats {
      flex-direction: column;
    }

    .stat {
      min-width: auto;
    }
  }
</style>
