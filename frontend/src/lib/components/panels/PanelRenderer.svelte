<script lang="ts">
  // Import simulation panels
  import {
    SpectrumPlot,
    IonFragmentation,
    StructurePanel,
    ConsoleOutput,
    PeakTable,
    ExportCenter,
  } from './simulation/index';

  // Import info panels
  import {
    SystemOverviewPanel,
    WorkflowPanel,
    MLEnginePanel,
    FragmentationSciencePanel,
    VisualizationFeaturesPanel,
    AccuracyValidationPanel,
    AdvancedFeaturesPanel,
    TechnicalRequirementsPanel,
  } from './info/index';

  // Import about panels
  import {
    DeveloperHeroPanel,
    DeveloperJourneyPanel,
    ProjectMetricsPanel,
    TechStackBanner,
    FeaturesShowcasePanel,
    OpenSourceCommitmentPanel,
    ContactConnectRevampedPanel,
    AcknowledgmentsStreamlinedPanel,
  } from './about/index';

  // Import chat components
  import { Chat as ChatWithSpectrum } from '../chat/index';
  import { CompoundStatsPanel } from './chat/index';

  // Import icons
  import { DynamicIcon } from '$lib/components/icons';
  import { devLog } from '$lib/utils/performance';

  export let component: string;
  export let props: Record<string, any> = {};
  export let isCarousel = false;

  // Log when props change
  $: {
    devLog.log(`PanelRenderer for ${component} props:`, props);
  }

  // Components mapping - expanded as we migrate components
  const components: Record<string, any> = {
    // Simulation panels
    SpectrumPlot,
    IonFragmentation,
    StructurePanel,
    ConsoleOutput,
    PeakTable,
    ExportCenter,

    // Chat components (✅ migrated)
    ChatWithSpectrum,
    CompoundStatsPanel,

    // Info panels (✅ migrated)
    SystemOverviewPanel,
    WorkflowPanel,
    MLEnginePanel,
    FragmentationSciencePanel,
    VisualizationFeaturesPanel,
    AccuracyValidationPanel,
    AdvancedFeaturesPanel,
    TechnicalRequirementsPanel,

    // About panels
    DeveloperHeroPanel,
    DeveloperJourneyPanel,
    ProjectMetricsPanel,
    TechStackBanner,
    FeaturesShowcasePanel,
    OpenSourceCommitmentPanel,
    ContactConnectRevampedPanel,
    AcknowledgmentsStreamlinedPanel,
  };

  $: ComponentToRender = components[component];
  $: enhancedProps = {
    ...props,
    ...(isCarousel ? { isCarousel: true } : {}),
    // Pass compact prop if it exists
    ...(props.compact !== undefined ? { compact: props.compact } : {}),
  };

  // Log enhanced props
  $: {
    devLog.log(`PanelRenderer enhanced props for ${component}:`, enhancedProps);
  }
</script>

{#if ComponentToRender}
  <svelte:component this={ComponentToRender} {...enhancedProps} />
{:else}
  <!-- Placeholder content for now -->
  <div class="placeholder">
    <div class="placeholder-content">
      <div class="placeholder-icon">
        <DynamicIcon name="SpectrumIcon" size={36} color="var(--text-secondary)" />
      </div>
      <h3>{component}</h3>
      <p>Panel component will be migrated here</p>
      {#if Object.keys(props).length > 0}
        <details>
          <summary>Props Debug</summary>
          <pre>{JSON.stringify(props, null, 2)}</pre>
        </details>
      {/if}
    </div>
  </div>
{/if}

<style>
  .placeholder {
    width: 100%;
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: rgba(255, 255, 255, 0.8);
    border-radius: var(--enforce-pill);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
    overflow: hidden !important;
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
    color: var(--text-secondary);
  }

  .placeholder h3 {
    margin: 0 0 0.5rem;
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text-primary);
    text-transform: capitalize;
  }

  .placeholder p {
    font-size: 0.9rem;
    color: var(--text-secondary);
    margin: 0 0 1rem;
  }

  details {
    font-size: 0.8rem;
    color: var(--text-tertiary);
    max-width: 300px;
  }

  details summary {
    cursor: pointer;
    margin-bottom: 0.5rem;
  }

  pre {
    text-align: left;
    background: rgba(0, 0, 0, 0.05);
    padding: 0.5rem;
    border-radius: 4px;
    overflow: auto;
    max-height: 200px;
  }
</style>
