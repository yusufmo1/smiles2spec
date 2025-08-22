<script lang="ts">
  import { onMount } from 'svelte';
  import Header from '$lib/components/Header.svelte';
  import SubtitlePill from '$lib/components/SubtitlePill.svelte';
  import { PanelGrid, PanelCarousel, PanelOverlay } from '$lib/components/panels/layout/index';
  import { currentPage, currentCarouselMode, focusedPanel, aboutPanels } from '$lib/stores/index';

  // Reactive values
  $: pageKey = $currentPage;
  $: carouselMode = $currentCarouselMode;
  $: overlayVisible = $focusedPanel !== null;
  $: currentPanels = $aboutPanels;

  function handleToggle(val: boolean) {
    if (val) {
      focusedPanel.set('developer-hero');
    } else {
      focusedPanel.set(null);
    }
  }

  function handleOverlayClose() {
    focusedPanel.set(null);
  }

  onMount(() => {
    console.log('About page mounted');
    console.log('Current page:', pageKey);
    console.log('Carousel mode:', carouselMode);
  });
</script>

<svelte:head>
  <title>About - Mass Spectrum Predictor</title>
  <meta
    name="description"
    content="Learn about the development, technology, and vision behind the Mass Spectrum Predictor platform."
  />
</svelte:head>

<div class="about-page">
  <Header
    className="prosto-one-regular"
    title="ABOUT"
    subtitle=""
    showToggle={true}
    toggleValue={carouselMode}
    onToggle={handleToggle}
  />

  <SubtitlePill
    icon="UserIcon"
    text="Learn about the development, technology, and vision behind this platform"
  />

  <div class="about-content">
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
  .about-page {
    width: 100%;
    min-height: 100vh;
  }

  .about-content {
    width: 100%;
    min-height: 60vh;
  }
</style>
