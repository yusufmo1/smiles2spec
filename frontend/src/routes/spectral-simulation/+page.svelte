<script lang="ts">
  import { onMount, onDestroy, tick } from 'svelte';
  import Header from '$lib/components/Header.svelte';
  import BulkNavigator from '$lib/components/BulkNavigator.svelte';
  import SubtitlePill from '$lib/components/SubtitlePill.svelte';
  import { SmilesInput } from '$lib/components/smiles-input/index';
  import { PanelGrid, PanelOverlay } from '$lib/components/panels/layout/index';
  import {
    currentPage,
    currentCarouselMode,
    currentCarouselIndex,
    focusedPanel,
    simulationPanels,
    appState,
    currentPredictionData,
    consoleText,
    initializePlots,
  } from '$lib/stores/index';
  import { predictSpectrum } from '$lib/services/api';
  import { DynamicIcon } from '$lib/components/icons';
  import { stateBatcher, devLog } from '$lib/utils/performance';

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
  $: bulkList = $appState.bulkList;
  $: bulkIndex = $appState.bulkIndex;
  $: consoleOutput = $consoleText;
  $: isLoading = $appState.isLoading;
  $: appError = $appState.error;

  // Reactive values
  $: pageKey = $currentPage;
  $: carouselMode = $currentCarouselMode;
  $: carouselIndex = $currentCarouselIndex;
  $: overlayVisible = $focusedPanel !== null;
  $: currentPanels = $simulationPanels;

  async function handleToggle(val: boolean) {
    if (val) {
      // Activate carousel mode by focusing first panel
      focusedPanel.set('spectrum');
      await tick();
    } else {
      // Deactivate by clearing focus
      focusedPanel.set(null);
    }
  }

  async function handlePredict(event: any) {
    const { smiles } = event.detail;

    try {
      // Batch initial state updates in single frame
      stateBatcher.queue(() => {
        appState.setLoading(true);
        appState.setError(null);
        appState.addConsoleEntry({
          type: 'info',
          message: `Starting prediction for: ${smiles}`,
          details: { smiles, action: 'predict_start' },
        });
      });

      const result = await predictSpectrum(smiles);
      
      // Batch success state updates in single frame
      stateBatcher.queue(() => {
        appState.setPredictionData(result);
        appState.addConsoleEntry({
          type: 'success',
          message: `Prediction completed successfully`,
          details: {
            smiles: result.smiles,
            peakCount: result.peaks?.length || 0,
            chemicalName: result.chemical_name,
            molecularWeight: result.molecular_weight,
            action: 'predict_success',
          },
        });
      });
    } catch (err: any) {
      // Batch error state updates in single frame
      stateBatcher.queue(() => {
        appState.setError(err.message);
        appState.addConsoleEntry({
          type: 'error',
          message: `Prediction failed: ${err.message}`,
          details: { error: err.message, smiles, action: 'predict_error' },
        });
      });
    } finally {
      // Batch final state updates
      stateBatcher.queue(() => {
        appState.setLoading(false);
        if (smilesInputRef) {
          smilesInputRef.setLoading(false);
        }
      });
    }
  }

  function handleBulk(event: any) {
    const { list } = event.detail;
    appState.setBulkList(list);
    appState.addConsoleEntry({
      type: 'info',
      message: `Bulk processing started: ${list.length} compounds`,
      details: { count: list.length, compounds: list, action: 'bulk_start' },
    });

    if (list.length > 0) {
      handlePredict({ detail: { smiles: list[0] } });
    }
  }

  function handleBulkNav(direction: 'prev' | 'next') {
    const newIndex = direction === 'next' ? $appState.bulkIndex + 1 : $appState.bulkIndex - 1;

    appState.setBulkIndex(newIndex);
    const smiles = $appState.bulkList[newIndex];

    appState.addConsoleEntry({
      type: 'info',
      message: `Navigated to compound ${newIndex + 1}/${$appState.bulkList.length}`,
      details: { index: newIndex, smiles, action: 'bulk_navigate' },
    });

    handlePredict({ detail: { smiles } });
  }

  function handleOverlayClose() {
    focusedPanel.set(null);
  }

  let plotCleanup: (() => void) | undefined;

  onMount(async () => {
    devLog.log('Spectral simulation page mounted');
    devLog.log('Current page:', pageKey);
    devLog.log('Carousel mode:', carouselMode);
    devLog.log('Carousel index:', carouselIndex);

    // Initialize plot system with lazy loading for better performance
    plotCleanup = await initializePlots();
  });

  onDestroy(() => {
    // Clean up plot effects when component is destroyed
    plotCleanup?.();
  });
</script>

<svelte:head>
  <title>Spectral Simulation - SMILES2SPEC</title>
  <meta
    name="description"
    content="Preview our foundational spectral simulation model. Predict high-resolution LC-ESI-qTOF mass spectra from SMILES notation using advanced machine learning."
  />
</svelte:head>

<div class="simulation-page">
  <Header
    className="prosto-one-regular"
    title="SPECTRAL SIMULATION"
    subtitle=""
    showToggle={hasFirstPrediction}
    toggleValue={carouselMode}
    onToggle={handleToggle}
  />

  <SubtitlePill
    icon="BeakerIcon"
    text="Preview our foundational model: Enter SMILES string to predict high-resolution LC-ESI-qTOF mass spectra"
  />

  <div class="simulation-content">
    <SmilesInput bind:this={smilesInputRef} on:predict={handlePredict} on:bulk={handleBulk} />

    {#if bulkList.length > 0}
      <BulkNavigator
        name={currentName}
        index={bulkIndex + 1}
        total={bulkList.length}
        on:prev={() => handleBulkNav('prev')}
        on:next={() => handleBulkNav('next')}
      />
    {/if}

    {#if appError}
      <div class="error-message">
        <span class="error-icon">
          <DynamicIcon name="ErrorIcon" size={20} color="#ff453a" />
        </span>
        {appError}
      </div>
    {/if}

    <div class="results-section">
      <PanelGrid
        currentPage={pageKey}
        {spectrumData}
        {peakData}
        {structurePNG}
        {currentSmiles}
        {currentName}
        consoleText={consoleOutput}
        {hasFirstPrediction}
        smilesList={bulkList}
      />
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
  .simulation-page {
    width: 100%;
    min-height: 100vh;
  }

  .simulation-content {
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

  .results-section {
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
    .simulation-content {
      gap: 1.5rem;
    }
  }
</style>
