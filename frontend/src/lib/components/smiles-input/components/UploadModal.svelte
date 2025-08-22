<script>
  import { createEventDispatcher } from 'svelte';
  import { smilesInputStore } from '../stores/smilesInputStore';
  import { smilesInputService } from '../services/smilesInputService';
  import DragDropZone from './DragDropZone.svelte';

  const dispatch = createEventDispatcher();

  let fileChooser;

  function closeModal() {
    smilesInputStore.hideUploadModal();
  }

  function handleEscapeKey(e) {
    if (e.key === 'Escape') {
      closeModal();
    }
  }

  async function handleFile(file) {
    if (!file) return;

    try {
      smilesInputStore.setLoading(true);
      const result = await smilesInputService.processFile(file);

      // Use store to trigger bulk upload instead of dispatching event
      smilesInputStore.triggerBulkUpload(result.smiles);
      closeModal();
    } catch (err) {
      alert(err.message);
    } finally {
      smilesInputStore.setLoading(false);
    }
  }

  function handleFileSelect(e) {
    const file = e.target.files[0];
    handleFile(file);
  }

  function handleSelectClick() {
    fileChooser.click();
  }
</script>

{#if $smilesInputStore.showUploadModal}
  <div
    class="modal-overlay"
    on:click={closeModal}
    on:keydown={handleEscapeKey}
    role="dialog"
    aria-modal="true"
    tabindex="0"
  >
    <div
      class="upload-modal glass-card"
      role="button"
      tabindex="0"
      on:click|stopPropagation
      on:keydown|stopPropagation
    >
      <button class="close-button" on:click={closeModal}>×</button>

      <div class="upload-content">
        <h3>Upload SMILES Data</h3>

        <DragDropZone on:file={(e) => handleFile(e.detail)}>
          <button slot="select-button" class="select-button" on:click={handleSelectClick}>
            Select File
          </button>
        </DragDropZone>
      </div>
    </div>
  </div>

  <input
    type="file"
    accept=".txt,.csv"
    bind:this={fileChooser}
    on:change={handleFileSelect}
    style="display:none"
  />
{/if}

<style>
  /* Modal Styles */
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

  .upload-modal {
    background: rgba(255, 255, 255, 0.9);
    width: 90%;
    max-width: 650px;
    min-height: 400px;
    padding: 2.5rem;
    border-radius: var(--enforce-pill);
    position: relative;
    box-shadow:
      0 15px 50px rgba(0, 0, 0, 0.2),
      0 0 0 1px rgba(255, 255, 255, 0.1);
    animation: slideUp 0.3s ease-out;
    overflow: hidden;
  }

  .upload-modal::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 30%;
    background: linear-gradient(to bottom, rgba(255, 255, 255, 0.5), rgba(255, 255, 255, 0.1));
    z-index: -1;
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

  .upload-content {
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
  }

  .upload-content h3 {
    margin: 0 0 1.5rem;
    font-size: 1.5rem;
    font-weight: 600;
    color: var(--text-primary);
  }

  .select-button {
    background: var(--accent);
    color: white;
    border: none;
    border-radius: var(--enforce-pill);
    padding: 0.75rem 1.5rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s;
    box-shadow: 0 4px 15px rgba(120, 121, 255, 0.3);
  }

  .select-button:hover {
    background: color-mix(in srgb, var(--accent), black 10%);
    transform: translateY(-2px);
  }

  .select-button:active {
    transform: translateY(0);
    box-shadow: 0 2px 8px rgba(120, 121, 255, 0.2);
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

  @media (max-width: 640px) {
    .modal-overlay {
      padding-top: 80px;
    }

    .upload-modal {
      width: 95%;
      max-width: none;
      min-height: 350px;
      margin: 0 1rem;
      padding: 1.5rem;
    }
  }
</style>
