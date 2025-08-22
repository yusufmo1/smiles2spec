<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
</script>

<div class="placeholder-state">
  <div class="placeholder-content">
    <div class="placeholder-icon">
      <DynamicIcon name="ChatWithSpectrumIcon" size={36} color="var(--text-secondary)" />
    </div>
    <h3>Chat with Spectrum</h3>
    <p>Predict a spectrum first to start chatting about molecular analysis</p>
  </div>
</div>

<style>
  .placeholder-state {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: rgba(255, 255, 255, 0.8);
    border-radius: var(--enforce-pill);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
  }

  .placeholder-content {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    opacity: 0.6;
    text-align: center;
    padding: 2rem;
  }

  .placeholder-icon {
    margin-bottom: 1rem;
    opacity: 0.7;
    display: flex;
  }

  .placeholder-content h3 {
    margin: 0 0 0.5rem;
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
  }

  .placeholder-content p {
    font-size: 0.9rem;
    color: var(--text-secondary);
    margin: 0;
  }
</style>
