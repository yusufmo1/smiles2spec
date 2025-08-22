<script lang="ts">
  import { onMount, onDestroy, tick } from 'svelte';
  import Header from '$lib/components/Header.svelte';
  import SubtitlePill from '$lib/components/SubtitlePill.svelte';
  import { SmilesInput } from '$lib/components/smiles-input/index';
  import { PanelGrid, PanelOverlay } from '$lib/components/panels/layout/index';
  import {
    currentPage,
    currentCarouselMode,
    focusedPanel,
    chatPanels,
    panelStore,
    appState,
    currentPredictionData,
    consoleText,
    initializePlots,
  } from '$lib/stores/index';
  import { predictSpectrum } from '$lib/services/api';
  import { DynamicIcon } from '$lib/components/icons';

  // State management - using global state
  let smilesInputRef: any;

  // Use global state instead of local state
  $: ({
    smiles: currentSmiles,
    name: currentName,
    spectrum: spectrumData,
    peaks: peakData,
    structure: structurePNG,
    molecularWeight,
    exactMass,
  } = $currentPredictionData);

  $: hasFirstPrediction = $appState.hasFirstPrediction;
  $: predictionConfidence = spectrumData ? 0.85 : 0; // Prediction confidence placeholder
  $: isLoading = $appState.isLoading;
  $: appError = $appState.error;

  // Reactive values
  $: pageKey = $currentPage;
  $: carouselMode = $currentCarouselMode;
  $: overlayVisible = $focusedPanel !== null;
  $: currentPanels = $chatPanels;

  // Chat panel navigation order
  const CHAT_PANEL_ORDER = ['chat-conversation', 'chat-spectrum', 'compound-stats'];

  // Automatically update panel props when global state changes
  // This will run whenever the prediction data changes, regardless of which page initiated it
  $: if ($currentPredictionData.smiles || hasFirstPrediction) {
    updatePanelProps();
  }

  // ✅ REPLACE: Simplified toggle handler
  async function handleToggle(val: boolean) {
    if (val) {
      // Activate carousel mode by focusing first panel
      focusedPanel.set('chat-conversation');
      await tick();
    } else {
      // Deactivate by clearing focus
      focusedPanel.set(null);
    }
  }

  async function handlePredict(event: any) {
    const { smiles } = event.detail;

    try {
      appState.setLoading(true);
      appState.setError(null);

      appState.addConsoleEntry({
        type: 'info',
        message: `Chat prediction started for: ${smiles}`,
        details: { smiles, source: 'chat_page', action: 'predict_start' },
      });

      const result = await predictSpectrum(smiles);
      appState.setPredictionData(result);

      appState.addConsoleEntry({
        type: 'success',
        message: `Chat prediction completed - Ready for analysis`,
        details: {
          smiles: result.smiles,
          peakCount: result.peaks?.length || 0,
          source: 'chat_page',
          action: 'predict_success',
        },
      });
    } catch (err: any) {
      appState.setError(err.message);
      appState.addConsoleEntry({
        type: 'error',
        message: `Chat prediction failed: ${err.message}`,
        details: { error: err.message, smiles, source: 'chat_page', action: 'predict_error' },
      });
    } finally {
      appState.setLoading(false);
      if (smilesInputRef) {
        smilesInputRef.setLoading(false);
      }
    }
  }

  function updatePanelProps() {
    // Update chat panel
    panelStore.updatePanelProps('chat-conversation', {
      hasSmilesPrediction: hasFirstPrediction,
      currentSmiles,
      isCarousel: false,
    });

    // Update spectrum panel
    panelStore.updatePanelProps('chat-spectrum', {
      spectrumData,
      isCarousel: false,
    });

    // Update compound stats panel
    panelStore.updatePanelProps('compound-stats', {
      smiles: currentSmiles,
      molecularWeight,
      exactMass,
      peakCount: peakData.length,
      chemicalName: currentName,
      predictionConfidence,
      isCarousel: false,
    });
  }

  function handleOverlayClose() {
    focusedPanel.set(null);
  }

  let plotCleanup: (() => void) | undefined;

  onMount(async () => {
    console.log('Chat with Spectrum page mounted');
    console.log('Current page:', pageKey);
    console.log('Carousel mode:', carouselMode);

    // Initialize plot system with lazy loading for better performance
    plotCleanup = await initializePlots();
  });

  onDestroy(() => {
    // Clean up plot effects when component is destroyed
    plotCleanup?.();
  });
</script>

<svelte:head>
  <title>Chat with Spectrum - Mass Spectrum Predictor</title>
  <meta
    name="description"
    content="Interactive chat interface for analysing and discussing mass spectra with AI assistance."
  />
</svelte:head>

<div class="chat-page">
  <Header
    className="prosto-one-regular"
    title="CHAT WITH SPECTRUM"
    subtitle=""
    showToggle={currentPanels.length > 1}
    toggleValue={carouselMode}
    onToggle={handleToggle}
  />

  <SubtitlePill
    icon="ChatWithSpectrumIcon"
    text="Interactive analysis and discussion of your mass spectra"
  />

  <div class="chat-content">
    <SmilesInput bind:this={smilesInputRef} on:predict={handlePredict} />

    {#if appError}
      <div class="error-message">
        <span class="error-icon">
          <DynamicIcon name="ErrorIcon" size={20} color="#ff453a" />
        </span>
        {appError}
      </div>
    {/if}

    <div class="chat-panels-section">
      <PanelGrid panels={currentPanels} currentPage={pageKey} />
    </div>
  </div>
</div>

<!-- Fullscreen Overlay -->
<PanelOverlay currentPage={pageKey} on:close={handleOverlayClose} />

<!-- Loading Overlay -->
{#if isLoading}
  <div class="loading-overlay">
    <div class="loading-content">
      <div class="loading-spinner">
        <DynamicIcon name="LoadingIcon" size={50} color="var(--accent)" />
      </div>
      <p>Analyzing spectrum...</p>
    </div>
  </div>
{/if}

<style>
  .chat-page {
    width: 100%;
    min-height: 100vh;
  }

  .chat-content {
    display: flex;
    flex-direction: column;
    gap: 2rem;
    padding-bottom: 2rem;
  }

  .error-message {
    background: rgba(255, 75, 85, 0.1);
    color: #ff453a;
    padding: 1rem 1.5rem;
    font-size: 0.9rem;
    font-weight: 500;
    display: flex;
    align-items: center;
    gap: 0.75rem;
    border-radius: 50px;
    border: 1px solid rgba(255, 75, 85, 0.25);
    box-shadow: 0 4px 12px rgba(255, 75, 85, 0.1);
    backdrop-filter: blur(10px);
  }

  .error-icon {
    flex-shrink: 0;
    display: flex;
    align-items: center;
    justify-content: center;
  }

  .chat-panels-section {
    width: 100%;
    min-height: 60vh;
  }

  .loading-overlay {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: rgba(247, 247, 248, 0.8);
    backdrop-filter: blur(10px);
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 9999;
  }

  .loading-content {
    text-align: center;
    color: var(--text-primary);
  }

  .loading-spinner {
    margin: 0 auto 1rem;
    display: flex;
    justify-content: center;
  }

  @media (max-width: 768px) {
    .chat-content {
      gap: 1.5rem;
    }
  }
</style>
