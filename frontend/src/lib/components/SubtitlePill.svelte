<script lang="ts">
  import { DynamicIcon, type IconName } from '$lib/components/icons';

  export let icon: IconName = 'BeakerIcon';
  export let text: string = '';
</script>

<div class="subtitle-pill">
  <DynamicIcon name={icon} size={16} color="var(--accent)" />
  <span>{text}</span>
</div>

<style>
  .subtitle-pill {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.5rem;
    padding: 0.5rem 1rem;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 20px;
    backdrop-filter: blur(10px);
    -webkit-backdrop-filter: blur(10px);
    border: 1px solid rgba(255, 255, 255, 0.1);
    font-size: 0.875rem;
    font-weight: 500;
    color: var(--text-primary);
    font-family:
      'SF Pro Display',
      system-ui,
      -apple-system,
      BlinkMacSystemFont,
      'Segoe UI',
      Roboto,
      sans-serif;
    transition: all var(--transition-smooth);
    margin: 0 auto 2rem;
    width: fit-content;
    margin-top: -1rem;
  }

  .subtitle-pill:hover {
    background: rgba(255, 255, 255, 0.12);
    transform: translateY(-1px);
  }
</style>
