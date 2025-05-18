<script>
  import { createEventDispatcher, onMount } from 'svelte';
  import { uploadSmilesFile } from '../services/api.js';
  
  const dispatch = createEventDispatcher();
  let smiles = '';
  let isLoading = false;
  let inputElement;
  let isFocused = false;
  let fileChooser;
  let showUploadModal = false;
  let dragActive = false;
  
  function handleSubmit() {
    if (!smiles.trim() || isLoading) return;
    isLoading = true;
    dispatch('predict', { smiles });
  }
  
  function handleKeydown(e) {
    if (e.key === 'Enter' && smiles.trim() && !isLoading) {
      handleSubmit();
    }
  }
  
  export function setLoading(loading) {
    isLoading = loading;
  }
  
  function openUploadModal() {
    showUploadModal = true;
  }
  
  function closeUploadModal() {
    showUploadModal = false;
    dragActive = false;
  }
  
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
      processFile(e.dataTransfer.files[0]);
    }
  }
  
  async function processFile(file) {
    const validTypes = ['text/plain', 'text/csv', 'application/csv', 'application/vnd.ms-excel'];
    const fileExt = file.name.split('.').pop().toLowerCase();
    
    if (!file || !(validTypes.includes(file.type) || ['txt', 'csv'].includes(fileExt))) {
      alert('Please upload a text (.txt) or CSV (.csv) file');
      return;
    }
    
    try {
      isLoading = true;
      const { smiles } = await uploadSmilesFile(file);
      dispatch('bulk', { list: smiles });
      closeUploadModal();
    } catch (err) {
      alert(err.message);
    } finally {
      isLoading = false;
    }
  }
  
  async function handleFile(e) {
    const file = e.target.files[0];
    if (!file) return;
    processFile(file);
  }
  
  onMount(() => {
    if (inputElement) inputElement.focus();
  });
</script>

<div class="input-container glass-card" class:focused={isFocused}>
  <div class="input-field">
    <span class="prompt-sign">›</span>
    <input 
      type="text" 
      bind:value={smiles}
      bind:this={inputElement}
      placeholder="Enter SMILES notation"
      disabled={isLoading}
      on:keydown={handleKeydown}
      on:focus={() => isFocused = true}
      on:blur={() => isFocused = false}
      aria-label="SMILES notation input"
    />
  </div>
  
  <!-- Upload button -->
  <button class="pill-button upload"
          aria-label="Upload SMILES list"
          on:click={openUploadModal}>
    <svg width="18" height="18" viewBox="0 0 24 24" fill="white">
      <path d="M12 3l4 4h-3v6h-2V7H8l4-4zm-7 8v9h14v-9h2v9a2 
               2 0 0 1-2 2H5a2 2 0 0 1-2-2v-9h2z"/>
    </svg>
  </button>
  
  <button 
    class="pill-button"
    on:click={handleSubmit} 
    disabled={!smiles.trim() || isLoading}
    aria-label="Analyze molecule"
  >
    {#if isLoading}
      <span class="loading-indicator"></span>
      <span>Processing</span>
    {:else}
      <span>Analyze</span>
    {/if}
  </button>
  
  <input type="file"
         accept=".txt,.csv"
         bind:this={fileChooser}
         on:change={handleFile}
         style="display:none" />
</div>

<!-- Upload Modal -->
{#if showUploadModal}
  <div class="modal-overlay" 
       on:click={closeUploadModal}
       on:keydown={(e) => e.key === 'Escape' && closeUploadModal()}
       role="dialog"
       aria-modal="true">
    <div class="upload-modal glass-card" 
         class:drag-active={dragActive}
         on:click|stopPropagation={() => {}}
         on:keydown|stopPropagation={() => {}}
         on:dragenter={handleDragEnter}
         on:dragover={handleDragOver}
         on:dragleave={handleDragLeave}
         on:drop={handleDrop}>
      
      <button class="close-button" on:click={closeUploadModal}>×</button>
      
      <div class="upload-content">
        <h3>Upload SMILES Data</h3>
        
        <div class="upload-area">
          <svg width="48" height="48" viewBox="0 0 24 24" fill="var(--accent)">
            <path d="M12 3l4 4h-3v6h-2V7H8l4-4zm-7 8v9h14v-9h2v9a2 
                   2 0 0 1-2 2H5a2 2 0 0 1-2-2v-9h2z"/>
          </svg>
          
          <p>Drag & drop a .txt or .csv file here<br> or</p>
          
          <button class="select-button" on:click={() => fileChooser.click()}>
            Select File
          </button>
          
          <p class="help-text">
            For TXT files: Each line should contain a single SMILES string.<br>
            For CSV files: First column should contain SMILES strings.
          </p>
        </div>
      </div>
    </div>
  </div>
{/if}

<style>
  .input-container {
    display: flex;
    align-items: center;
    margin: 0 0 2rem;
    padding: 0.65rem;
    transition: all var(--transition-elastic);
    background: var(--surface-glass-hover);
    border-radius: var(--enforce-pill);
    position: relative;
    overflow: hidden !important;
  }
  
  .input-container.focused {
    box-shadow: 
      0 10px 35px 0 rgba(0, 0, 0, 0.08),
      0 0 0 2px var(--accent-soft),
      0 2px 8px 0 rgba(255, 255, 255, 0.4) inset;
  }
  
  /* Highlight effect for focused state */
  .input-container.focused::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 40%;
    background: linear-gradient(
      to bottom,
      rgba(255, 255, 255, 0.25),
      rgba(255, 255, 255, 0.05)
    );
    z-index: -1;
  }
  
  .input-field {
    display: flex;
    align-items: center;
    flex: 1;
    padding-left: 1.25rem;
  }
  
  .prompt-sign {
    color: var(--accent);
    margin-right: 0.75rem;
    font-weight: 500;
    font-size: 1.25rem;
    opacity: 0.85;
    text-shadow: 0 0 8px rgba(120, 121, 255, 0.3);
  }
  
  input {
    flex: 1;
    background: transparent;
    border: none;
    color: var(--text-primary);
    font-family: inherit;
    font-size: 1rem;
    padding: 0.85rem 0;
    outline: none;
  }
  
  input::placeholder {
    color: var(--text-tertiary);
  }
  
  input:disabled {
    opacity: 0.7;
    cursor: not-allowed;
  }
  
  .pill-button {
    margin-right: 0.25rem;
    padding: 0.85rem 1.8rem;
    border-radius: var(--enforce-pill);
    transition: all var(--transition-elastic);
    font-weight: 500;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    color: white;
    border: none;
    font-size: 0.95rem;
    cursor: pointer;
    box-shadow: 
      0 6px 18px rgba(120, 121, 255, 0.25),
      0 0 0 1px rgba(255, 255, 255, 0.1);
  }
  
  .upload {
    background: var(--accent-secondary);
    padding: 0.65rem;
    display: flex;
    align-items: center;
    justify-content: center;
  }
  
  .upload:hover { background: var(--accent); }
  .upload svg { pointer-events: none; }
  
  .pill-button:hover:not(:disabled) {
    background: linear-gradient(135deg, var(--accent-secondary) 0%, var(--accent) 100%);
    transform: translateY(-2px) scale(1.025);
    box-shadow: 
      0 8px 24px rgba(120, 121, 255, 0.35),
      0 0 0 1px rgba(255, 255, 255, 0.15);
  }
  
  .pill-button:disabled {
    background: rgba(0, 0, 0, 0.1);
    color: var(--text-tertiary);
    cursor: not-allowed;
    box-shadow: none;
  }
  
  /* Loading animation */
  .loading-indicator {
    width: 1.1rem;
    height: 1.1rem;
    border: 2px solid rgba(255, 255, 255, 0.3);
    border-top: 2px solid white;
    border-radius: 50%;
    animation: spin 1s linear infinite;
  }
  
  @keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
  }
  
  /* Modal Styles */
  .modal-overlay {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: rgba(0, 0, 0, 0.5);
    backdrop-filter: blur(5px);
    display: flex;
    justify-content: center;
    align-items: flex-start;
    padding-top: 120px;
    z-index: 1000;
    animation: fadeIn 0.2s ease-out;
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
    background: linear-gradient(
      to bottom,
      rgba(255, 255, 255, 0.5),
      rgba(255, 255, 255, 0.1)
    );
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
  }
  
  .drag-active .upload-area {
    border-color: var(--accent);
    background: rgba(120, 121, 255, 0.05);
    transform: scale(1.02);
  }
  
  .upload-area p {
    margin: 1rem 0;
    color: var(--text-secondary);
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
  
  @keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
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
    .input-container {
      flex-direction: column;
      padding: 1.25rem;
    }
    
    .input-field {
      width: 100%;
      margin-bottom: 1.25rem;
      padding-left: 0;
    }
    
    .pill-button {
      width: 100%;
      margin-left: 0;
      justify-content: center;
    }
    
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
    
    .upload-area {
      padding: 2rem 1rem;
      min-height: 200px;
    }
  }
</style> 