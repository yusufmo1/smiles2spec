<script>
  import SpectrumPlot from './SpectrumPlot.svelte';
  import IonFragmentation from './IonFragmentation.svelte';
  import StructurePanel from './StructurePanel.svelte';
  import PeakTable from './PeakTable.svelte';
  import ConsoleOutput from './ConsoleOutput.svelte';
  import ChatWithSpectrum from './ChatWithSpectrum.svelte';

  export let component;
  export let props = {};
  export let isCarousel = false;

  const components = {
    SpectrumPlot,
    IonFragmentation,
    StructurePanel,
    PeakTable,
    ConsoleOutput,
    ChatWithSpectrum
  };

  $: ComponentToRender = components[component];
  $: enhancedProps = isCarousel ? { ...props, isCarousel: true } : props;
</script>

{#if ComponentToRender}
  <svelte:component this={ComponentToRender} {...enhancedProps} />
{:else}
  <div class="error">Component '{component}' not found</div>
{/if}

<style>
  .error {
    display: flex;
    align-items: center;
    justify-content: center;
    height: 100%;
    color: var(--text-secondary);
    font-style: italic;
  }
</style>
