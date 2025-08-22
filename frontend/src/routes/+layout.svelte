<script lang="ts">
  import '../app.css';
  import { page } from '$app/stores';
  import Background from '$lib/components/Background.svelte';
  import ErrorBoundary from '$lib/components/ErrorBoundary.svelte';
  import Navbar from '$lib/components/Navbar.svelte';
  import { focusedPanel } from '$lib/stores/index';

  $: currentPath = $page.url.pathname;
  $: isFullScreenCarousel = $focusedPanel !== null;

  function handleApplicationError(event: any) {
    console.error('🚨 Application Error Captured:', event.detail);
    // Could send error to logging service here
  }
</script>

<svelte:head>
  <link rel="preconnect" href="https://fonts.googleapis.com" />
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin="anonymous" />
  <link
    href="https://fonts.googleapis.com/css2?family=SF+Pro+Display:wght@400;500;600&family=Prosto+One:wght@400&display=swap"
    rel="stylesheet"
  />
</svelte:head>

<!-- Add the animated background -->
<Background />

<!-- Add error boundary -->
<ErrorBoundary on:error={handleApplicationError}>
  <!-- Navigation - HIDE IN CAROUSEL MODE -->
  {#if !isFullScreenCarousel}
    <Navbar />
  {/if}

  <!-- App shell with conditional styling -->
  <div class="app-shell" class:fullscreen={isFullScreenCarousel}>
    <slot />
  </div>
</ErrorBoundary>

<!-- Note: PanelOverlay is handled in individual page components -->
<!-- Note: Modals are handled by individual page components -->

<style>
  /* App shell styling */
  .app-shell {
    width: 100%;
    max-width: 90%;
    margin: 2rem auto 3rem;
    padding: 2rem;
    border-radius: 24px;
    background: rgba(255, 255, 255, 0.35);
    box-shadow:
      0 25px 55px rgba(0, 0, 0, 0.08),
      0 3px 12px rgba(0, 0, 0, 0.04);
    backdrop-filter: blur(30px) saturate(160%);
    -webkit-backdrop-filter: blur(30px) saturate(160%);
    min-height: 70vh;
    transition: all 0.3s ease;

    /* Ensure this doesn't create stacking context issues */
    position: relative;
    z-index: 1;
  }

  /* FULL-SCREEN CAROUSEL MODE */
  .app-shell.fullscreen {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    max-width: 100%;
    margin: 0;
    padding: 0;
    border-radius: 0;
    background: rgba(247, 247, 248, 0.95);
    backdrop-filter: blur(25px) saturate(180%);
    z-index: 9999;
    min-height: 100vh;
    height: 100vh;
    overflow: hidden;
  }

  /* Ensure fullscreen overlays are not constrained */
  :global(.fullscreen-overlay) {
    position: fixed !important;
    top: 0 !important;
    left: 0 !important;
    right: 0 !important;
    bottom: 0 !important;
    z-index: 10000 !important;
  }

  /* Responsive adjustments */
  @media (max-width: 768px) {
    .app-shell:not(.fullscreen) {
      max-width: 95%;
      padding: 1.5rem;
    }
  }
</style>
