<script>
  import { smilesInputStore } from '../stores/smilesInputStore';
  import { smilesInputService } from '../services/smilesInputService';
  import { DynamicIcon } from '$lib/components/icons';

  let smilesCount = 1;
  let smilesDescription = '';
  let isGenerating = false;

  function closeModal() {
    smilesInputStore.hideGenerateModal();
  }

  function handleEscapeKey(e) {
    if (e.key === 'Escape') {
      closeModal();
    }
  }

  async function generateRandomSmiles() {
    try {
      isGenerating = true;
      const result = await smilesInputService.generateSmiles({
        count: smilesCount,
        description: smilesDescription,
      });

      // For multiple SMILES, join with newlines
      if (result.smiles.length > 0) {
        smilesInputStore.setSmiles(result.smiles.join('\n'));
      }

      closeModal();

      // No auto-submission - user must manually click Analyze
    } catch (err) {
      alert(err.message);
    } finally {
      isGenerating = false;
    }
  }
</script>

{#if $smilesInputStore.showGenerateModal}
  <div
    class="modal-overlay"
    on:click={closeModal}
    on:keydown={handleEscapeKey}
    role="dialog"
    tabindex="-1"
    aria-modal="true"
  >
    <div
      class="generate-modal glass-card"
      role="button"
      tabindex="0"
      on:click|stopPropagation
      on:keydown|stopPropagation
    >
      <button class="close-button" on:click={closeModal}>×</button>

      <div class="generate-content">
        <h3>Generate Random SMILES</h3>

        <div class="generate-area">
          <div class="wand-icon-large">
            <DynamicIcon name="WandIcon" size={48} color="#9c27b0" />
          </div>

          <p>Generate a random molecule with valid SMILES notation</p>

          <div class="slider-container">
            <label for="smiles-count"
              >Number of SMILES to generate: <span class="count-value">{smilesCount}</span></label
            >
            <div class="custom-slider-container">
              <input
                type="range"
                id="smiles-count"
                min="1"
                max="10"
                step="1"
                bind:value={smilesCount}
                class="custom-slider"
              />
              <div class="slider-markers">
                {#each Array(10) as _, i}
                  <div class="marker-value">{i + 1}</div>
                {/each}
              </div>
            </div>
          </div>

          <div class="description-container">
            <label for="smiles-description">Describe the type of molecule (optional):</label>
            <input
              type="text"
              id="smiles-description"
              placeholder="e.g., 'a drug-like molecule with a benzene ring'"
              bind:value={smilesDescription}
            />
          </div>

          <button class="generate-button" on:click={generateRandomSmiles} disabled={isGenerating}>
            {#if isGenerating}
              <span class="loading-indicator"></span>
              <span>Generating...</span>
            {:else}
              <span>Generate {smilesCount > 1 ? `${smilesCount} SMILES` : 'Random SMILES'}</span>
            {/if}
          </button>
        </div>
      </div>
    </div>
  </div>
{/if}

<style>
  .modal-overlay {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    width: 100vw;
    height: 100vh;
    background: transparent;
    backdrop-filter: blur(15px);
    display: flex;
    justify-content: center;
    align-items: center;
    padding-top: 0;
    z-index: 20000;
    animation: overlay-fade-in 0.2s ease-out;
  }

  .generate-modal {
    background: rgba(255, 255, 255, 0.9);
    width: 90%;
    max-width: 650px;
    min-height: 350px;
    padding: 2.5rem;
    border-radius: var(--enforce-pill);
    position: relative;
    box-shadow:
      0 15px 50px rgba(0, 0, 0, 0.2),
      0 0 0 1px rgba(255, 255, 255, 0.1);
    animation: slideUp 0.3s ease-out;
    overflow: hidden;
  }

  .close-button {
    position: absolute;
    top: 1rem;
    right: 1.5rem;
    font-size: 1.5rem;
    width: 2rem;
    height: 2rem;
    display: flex;
    align-items: center;
    justify-content: center;
    border: none;
    background: transparent;
    cursor: pointer;
    color: var(--text-secondary);
    transition: all 0.2s;
  }

  .close-button:hover {
    color: var(--text-primary);
    transform: scale(1.1);
  }

  .generate-content {
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
  }

  .generate-content h3 {
    margin: 0 0 1.5rem;
    font-size: 1.5rem;
    font-weight: 600;
    color: var(--text-primary);
  }

  .generate-area {
    border: 2px dashed #9c27b0;
    border-radius: 1rem;
    padding: 3rem;
    width: 100%;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 1.5rem;
  }

  .wand-icon-large {
    margin-bottom: 1rem;
  }

  .slider-container {
    width: 100%;
    max-width: 400px;
  }

  .custom-slider-container {
    margin-top: 1.5rem;
    position: relative;
    padding: 0 10px;
  }

  .custom-slider {
    -webkit-appearance: none;
    appearance: none;
    width: 100%;
    height: 8px;
    border-radius: 4px;
    background: linear-gradient(to right, #e1bee7, #9c27b0);
    outline: none;
    margin: 15px 0;
  }

  .custom-slider::-webkit-slider-thumb {
    -webkit-appearance: none;
    appearance: none;
    width: 25px;
    height: 25px;
    border-radius: 50%;
    background: #9c27b0;
    cursor: pointer;
    border: 3px solid white;
    box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2);
    transition: all 0.2s ease;
  }

  .custom-slider::-moz-range-thumb {
    width: 25px;
    height: 25px;
    border-radius: 50%;
    background: #9c27b0;
    cursor: pointer;
    border: 3px solid white;
    box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2);
    transition: all 0.2s ease;
  }

  .custom-slider::-webkit-slider-thumb:hover {
    background: #7b1fa2;
    transform: scale(1.2);
  }

  .custom-slider::-moz-range-thumb:hover {
    background: #7b1fa2;
    transform: scale(1.2);
  }

  .slider-markers {
    display: flex;
    justify-content: space-between;
    width: 100%;
    padding: 0 10px;
    margin-top: -8px;
  }

  .marker-value {
    font-size: 0.8rem;
    color: var(--text-secondary);
    position: relative;
    display: flex;
    justify-content: center;
    width: 20px;
    margin-left: -10px;
  }

  .count-value {
    font-weight: bold;
    color: #9c27b0;
  }

  .description-container {
    width: 100%;
    max-width: 400px;
  }

  .description-container input {
    width: 100%;
    padding: 0.75rem;
    border: 1px solid var(--surface-stroke);
    border-radius: var(--radius-md);
    margin-top: 0.5rem;
  }

  .generate-button {
    background: #9c27b0;
    color: white;
    border: none;
    border-radius: var(--enforce-pill);
    padding: 0.75rem 1.5rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s;
    display: flex;
    align-items: center;
    gap: 0.5rem;
  }

  .generate-button:hover:not(:disabled) {
    background: #7b1fa2;
    transform: translateY(-2px);
  }

  .generate-button:disabled {
    opacity: 0.6;
    cursor: not-allowed;
  }

  .loading-indicator {
    width: 1rem;
    height: 1rem;
    border: 2px solid rgba(255, 255, 255, 0.3);
    border-top: 2px solid white;
    border-radius: 50%;
    animation: spin 1s linear infinite;
  }

  @keyframes spin {
    0% {
      transform: rotate(0deg);
    }
    100% {
      transform: rotate(360deg);
    }
  }

  @keyframes overlay-fade-in {
    from {
      opacity: 0;
      backdrop-filter: blur(0px);
    }
    to {
      opacity: 1;
      backdrop-filter: blur(15px);
    }
  }

  @keyframes slideUp {
    from {
      opacity: 0;
      transform: translateY(30px);
    }
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }
</style>
