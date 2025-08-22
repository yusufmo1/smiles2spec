<script lang="ts">
  import { onMount } from 'svelte';
  import Header from '$lib/components/Header.svelte';
  import SubtitlePill from '$lib/components/SubtitlePill.svelte';
  import { PanelGrid, PanelCarousel, PanelOverlay } from '$lib/components/panels/layout/index';
  import { currentPage, currentCarouselMode, focusedPanel, infoPanels } from '$lib/stores/index';

  // Reactive values
  $: pageKey = $currentPage;
  $: carouselMode = $currentCarouselMode;
  $: overlayVisible = $focusedPanel !== null;
  $: currentPanels = $infoPanels;

  function handleToggle(val: boolean) {
    if (val) {
      focusedPanel.set('overview');
    } else {
      focusedPanel.set(null);
    }
  }

  function handleOverlayClose() {
    focusedPanel.set(null);
  }

  onMount(() => {
    console.log('How it works page mounted');
    console.log('Current page:', pageKey);
    console.log('Carousel mode:', carouselMode);
  });
</script>

<svelte:head>
  <title>How It Works - Mass Spectrum Predictor</title>
  <meta
    name="description"
    content="Understand the science and technology behind our foundational spectral simulation model and prediction methodologies."
  />
</svelte:head>

<div class="how-it-works-page">
  <Header
    className="prosto-one-regular"
    title="HOW IT WORKS"
    subtitle=""
    showToggle={true}
    toggleValue={carouselMode}
    onToggle={handleToggle}
  />

  <SubtitlePill
    icon="LightbulbIcon"
    text="Understand the science and technology behind our foundational spectral simulation model"
  />

  <div class="how-it-works-content">
    {#if carouselMode}
      <PanelCarousel {currentPanels} />
    {:else}
      <PanelGrid panels={currentPanels} currentPage={pageKey} />
    {/if}
  </div>
</div>

<!-- Fullscreen Overlay -->
<PanelOverlay currentPage={pageKey} on:close={handleOverlayClose} />

<style>
  .how-it-works-page {
    width: 100%;
    min-height: 100vh;
  }

  .how-it-works-content {
    width: 100%;
    min-height: 60vh;
  }
</style>
