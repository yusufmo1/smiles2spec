<script>
  import { createEventDispatcher } from 'svelte';

  const dispatch = createEventDispatcher();

  let dragActive = false;

  function handleDragEnter(e) {
    e.preventDefault();
    e.stopPropagation();
    dragActive = true;
  }

  function handleDragOver(e) {
    e.preventDefault();
    e.stopPropagation();
    if (!dragActive) dragActive = true;
  }

  function handleDragLeave(e) {
    e.preventDefault();
    e.stopPropagation();
    dragActive = false;
  }

  function handleDrop(e) {
    e.preventDefault();
    e.stopPropagation();
    dragActive = false;

    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      dispatch('file', e.dataTransfer.files[0]);
    }
  }
</script>

<div
  class="upload-area"
  class:drag-active={dragActive}
  on:dragenter={handleDragEnter}
  on:dragover={handleDragOver}
  on:dragleave={handleDragLeave}
  on:drop={handleDrop}
  role="button"
  tabindex="0"
  aria-label="Drop file here"
>
  <svg width="48" height="48" viewBox="0 0 24 24" fill="var(--accent)">
    <path
      d="M12 3l4 4h-3v6h-2V7H8l4-4zm-7 8v9h14v-9h2v9a2 
             2 0 0 1-2 2H5a2 2 0 0 1-2-2v-9h2z"
    />
  </svg>

  <p>Drag & drop a .txt or .csv file here<br /> or</p>

  <slot name="select-button">
    <button class="select-button" on:click> Select File </button>
  </slot>

  <p class="help-text">
    For TXT files: Each line should contain a single SMILES string.<br />
    For CSV files: First column should contain SMILES strings.
  </p>
</div>

<style>
  .upload-area {
    border: 2px dashed var(--accent-soft);
    border-radius: 1rem;
    padding: 3rem;
    width: 100%;
    min-height: 240px;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    transition: all 0.2s;
    cursor: pointer;
  }

  .drag-active {
    border-color: var(--accent);
    background: rgba(120, 121, 255, 0.05);
    transform: scale(1.02);
  }

  .upload-area p {
    margin: 1rem 0;
    color: var(--text-secondary);
    text-align: center;
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
    background: var(--accent-secondary);
    transform: translateY(-2px);
    box-shadow: 0 6px 20px rgba(120, 121, 255, 0.4);
  }

  .help-text {
    font-size: 0.85rem;
    color: var(--text-tertiary);
    margin-top: 1.5rem;
  }

  @media (max-width: 640px) {
    .upload-area {
      padding: 2rem 1rem;
      min-height: 200px;
    }
  }
</style>
