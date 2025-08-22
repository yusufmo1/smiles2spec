<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';

  export let steps: Array<{
    id: string;
    title: string;
    icon: any;
    color: string;
    description?: string;
  }> = [];
  export let direction: 'horizontal' | 'vertical' = 'horizontal';

  let hoveredStep: string | null = null;
</script>

<div class="process-flow" class:vertical={direction === 'vertical'}>
  {#each steps as step, i}
    <div
      class="flow-step"
      class:hovered={hoveredStep === step.id}
      on:mouseenter={() => (hoveredStep = step.id)}
      on:mouseleave={() => (hoveredStep = null)}
      style="--step-color: {step.color}; --step-delay: {i * 0.1}s"
      role="button"
      tabindex="0"
      on:keydown={(e) =>
        e.key === 'Enter' && (hoveredStep = hoveredStep === step.id ? null : step.id)}
    >
      <div class="step-icon">
        <div class="icon-glow"></div>
        <DynamicIcon name={step.icon as any} size={32} color="white" />
      </div>
      <h4>{step.title}</h4>
      {#if step.description && hoveredStep === step.id}
        <p class="step-description">{step.description}</p>
      {/if}
    </div>

    {#if i < steps.length - 1}
      <div class="flow-connector">
        <div class="connector-line"></div>
        <div class="connector-pulse"></div>
      </div>
    {/if}
  {/each}
</div>

<style>
  .process-flow {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1rem;
    height: 100%;
    padding: 1rem;
  }

  .process-flow.vertical {
    flex-direction: column;
  }

  .flow-step {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
    position: relative;
    cursor: pointer;
    transition: all 0.3s ease;
    animation: fadeInUp 0.6s ease-out forwards;
    animation-delay: var(--step-delay);
    opacity: 0;
  }

  @keyframes fadeInUp {
    to {
      opacity: 1;
      transform: translateY(0);
    }
    from {
      opacity: 0;
      transform: translateY(20px);
    }
  }

  .flow-step:hover {
    transform: translateY(-5px);
  }

  .step-icon {
    width: 80px;
    height: 80px;
    background: var(--step-color);
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    margin-bottom: 1rem;
    position: relative;
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.15);
    transition: all 0.3s ease;
  }

  .flow-step:hover .step-icon {
    transform: scale(1.1);
    box-shadow: 0 12px 32px rgba(0, 0, 0, 0.2);
  }

  .icon-glow {
    position: absolute;
    inset: -20px;
    background: radial-gradient(circle, var(--step-color), transparent 70%);
    opacity: 0;
    transition: opacity 0.3s ease;
  }

  .flow-step:hover .icon-glow {
    opacity: 0.3;
  }

  h4 {
    margin: 0 0 0.5rem;
    font-size: 1rem;
    color: var(--text-primary);
    font-weight: 600;
  }

  .step-description {
    position: absolute;
    top: 100%;
    left: 50%;
    transform: translateX(-50%);
    background: rgba(0, 0, 0, 0.9);
    color: white;
    padding: 0.75rem 1rem;
    border-radius: 8px;
    font-size: 0.85rem;
    white-space: nowrap;
    margin-top: 0.5rem;
    opacity: 0;
    animation: fadeIn 0.2s ease-out forwards;
    z-index: 10;
  }

  @keyframes fadeIn {
    to {
      opacity: 1;
    }
  }

  .flow-connector {
    width: 60px;
    height: 2px;
    position: relative;
    align-self: center;
    margin: 0 -0.5rem;
  }

  .connector-line {
    width: 100%;
    height: 100%;
    background: linear-gradient(90deg, var(--accent), var(--accent-secondary));
    opacity: 0.3;
  }

  .connector-pulse {
    position: absolute;
    top: 50%;
    left: 0;
    width: 20px;
    height: 20px;
    background: var(--accent);
    border-radius: 50%;
    transform: translate(-50%, -50%);
    animation: pulse 2s linear infinite;
  }

  @keyframes pulse {
    0% {
      transform: translate(-50%, -50%) scale(0);
      opacity: 1;
    }
    100% {
      transform: translate(250%, -50%) scale(0);
      opacity: 0;
    }
  }

  .process-flow.vertical .flow-connector {
    width: 2px;
    height: 60px;
    margin: -0.5rem 0;
  }

  .process-flow.vertical .connector-pulse {
    animation: pulseVertical 2s linear infinite;
  }

  @keyframes pulseVertical {
    0% {
      transform: translate(-50%, -50%) scale(0);
      opacity: 1;
    }
    100% {
      transform: translate(-50%, 250%) scale(0);
      opacity: 0;
    }
  }
</style>
