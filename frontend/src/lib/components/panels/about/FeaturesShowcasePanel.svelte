<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let compact: boolean = false;
</script>

<div class="features-showcase" class:compact>
  <div class="bento-grid">
    <div class="feature-card large-card">
      <div class="feature-icon">
        <DynamicIcon name="SpectrumIcon" size={compact ? 32 : 48} />
      </div>
      <h3 class="feature-title">Real-time Visualization</h3>
      <p class="feature-description">
        Interactive spectrum visualization with zoom, pan, and fragment highlighting.
      </p>
      <div class="mini-spectrum">
        <div class="spectrum-line"></div>
        <div class="spectrum-line"></div>
        <div class="spectrum-line"></div>
        <div class="spectrum-line"></div>
        <div class="spectrum-line"></div>
      </div>
    </div>

    <div class="feature-card medium-card">
      <div class="feature-icon">
        <DynamicIcon name="ChatWithSpectrumIcon" size={compact ? 28 : 40} />
      </div>
      <h3 class="feature-title">AI Chat Integration</h3>
      <p class="feature-description">
        Chat with your spectrum data to get insights and interpretations.
      </p>
    </div>

    <div class="feature-card small-card">
      <div class="feature-icon">
        <DynamicIcon name="DownloadIcon" size={compact ? 24 : 32} />
      </div>
      <h3 class="feature-title">Export Options</h3>
      <p class="feature-description">Multiple formats for easy sharing and publication.</p>
    </div>

    <div class="feature-card small-card">
      <div class="feature-icon">
        <DynamicIcon name="FlashIcon" size={compact ? 24 : 32} />
      </div>
      <h3 class="feature-title">Bulk Processing</h3>
      <p class="feature-description">Process multiple structures simultaneously.</p>
    </div>

    <div class="feature-card medium-card">
      <div class="feature-icon">
        <DynamicIcon name="SpectrumIcon" size={compact ? 28 : 40} />
      </div>
      <h3 class="feature-title">Advanced Analytics</h3>
      <p class="feature-description">
        Deeper insights with peak analysis and structural correlation.
      </p>
    </div>
  </div>
</div>

<style>
  .features-showcase {
    display: flex;
    justify-content: center;
    padding: 3rem 2rem;
  }

  .features-showcase.compact {
    padding: 2rem 1rem;
  }

  .bento-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    grid-template-rows: repeat(3, auto);
    gap: 1.5rem;
    width: 100%;
    max-width: 1200px;
  }

  .feature-card {
    position: relative;
    background: var(--surface-glass);
    border-radius: var(--radius-lg);
    padding: 2rem;
    border: 1px solid rgba(255, 255, 255, 0.2);
    display: flex;
    flex-direction: column;
    transition:
      transform 0.3s ease,
      box-shadow 0.3s ease;
    overflow: hidden;
  }

  .feature-card:hover {
    transform: translateY(-5px);
    box-shadow: 0 15px 30px rgba(0, 0, 0, 0.1);
  }

  .large-card {
    grid-column: span 2;
    grid-row: span 2;
  }

  .medium-card {
    grid-column: span 1;
    grid-row: span 2;
  }

  .small-card {
    grid-column: span 1;
    grid-row: span 1;
  }

  .feature-icon {
    color: var(--accent);
    margin-bottom: 1.5rem;
    display: flex;
    align-items: center;
  }

  .feature-title {
    font-size: 1.2rem;
    font-weight: 600;
    margin: 0 0 0.75rem 0;
    color: var(--text-primary);
  }

  .feature-description {
    font-size: 0.95rem;
    line-height: 1.6;
    color: var(--text-secondary);
    margin: 0;
    flex-grow: 1;
  }

  .mini-spectrum {
    margin-top: 2rem;
    height: 100px;
    display: flex;
    align-items: flex-end;
    gap: 8px;
  }

  .spectrum-line {
    width: 6px;
    background: var(--accent);
    border-radius: 4px 4px 0 0;
    position: relative;
    animation: pulse 2s infinite alternate;
  }

  .spectrum-line:nth-child(1) {
    height: 40%;
    animation-delay: 0s;
  }

  .spectrum-line:nth-child(2) {
    height: 70%;
    animation-delay: 0.2s;
  }

  .spectrum-line:nth-child(3) {
    height: 100%;
    animation-delay: 0.4s;
  }

  .spectrum-line:nth-child(4) {
    height: 60%;
    animation-delay: 0.6s;
  }

  .spectrum-line:nth-child(5) {
    height: 30%;
    animation-delay: 0.8s;
  }

  @keyframes pulse {
    0% {
      opacity: 0.6;
    }
    100% {
      opacity: 1;
      box-shadow: 0 0 10px rgba(120, 121, 255, 0.8);
    }
  }

  @media (max-width: 992px) {
    .bento-grid {
      grid-template-columns: repeat(2, 1fr);
      grid-template-rows: auto;
    }

    .large-card {
      grid-column: span 2;
    }

    .medium-card {
      grid-column: span 1;
    }

    .small-card {
      grid-column: span 1;
    }
  }

  @media (max-width: 768px) {
    .bento-grid {
      grid-template-columns: 1fr;
    }

    .large-card,
    .medium-card,
    .small-card {
      grid-column: span 1;
      grid-row: auto;
    }

    .feature-card {
      padding: 1.5rem;
    }
  }
</style>
