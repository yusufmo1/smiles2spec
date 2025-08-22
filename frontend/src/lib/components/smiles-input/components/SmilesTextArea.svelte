<script>
  import { createEventDispatcher, onMount } from 'svelte';
  import { smilesInputStore, lineCount } from '../stores/smilesInputStore';

  const dispatch = createEventDispatcher();

  let textareaElement;
  let prefixEl;

  // Handle keyboard shortcuts
  function handleKeydown(e) {
    // Only submit on Ctrl+Enter or Cmd+Enter
    if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') {
      e.preventDefault();
      dispatch('submit');
    }
  }

  // Keep arrow prefix in sync with textarea content
  function syncPrefix() {
    if (!prefixEl) return;

    const rows = Math.max(1, $lineCount);
    prefixEl.textContent = Array(rows).fill('›').join('\n');
  }

  // Reactive statement to sync prefix when smiles or lineCount changes
  $: {
    $smilesInputStore.smiles; // Track smiles changes
    $lineCount; // Track line count changes
    if (prefixEl) syncPrefix();
  }

  // Focus management
  function handleFocus() {
    smilesInputStore.setFocused(true);
  }

  function handleBlur() {
    smilesInputStore.setFocused(false);
  }

  // Input handling
  function handleInput() {
    // Sync prefix on every input change
    syncPrefix();
  }

  onMount(() => {
    if (textareaElement) {
      textareaElement.focus();
    }
    // Initialize the arrow prefix
    syncPrefix();
  });
</script>

<div class="multi-input-wrapper">
  <pre class="line-prefix" bind:this={prefixEl}>›</pre>
  <textarea
    bind:value={$smilesInputStore.smiles}
    bind:this={textareaElement}
    placeholder="Enter SMILES notation"
    disabled={$smilesInputStore.isLoading}
    rows={Math.min(8, $lineCount)}
    on:input={handleInput}
    on:keydown={handleKeydown}
    on:focus={handleFocus}
    on:blur={handleBlur}
    aria-label="SMILES notation input"
  ></textarea>
</div>

<style>
  /* Multi-line input wrapper */
  .multi-input-wrapper {
    position: relative;
    flex: 1;
  }

  /* Arrow prefix styling */
  .line-prefix {
    position: absolute;
    top: 0;
    left: 0;
    padding: 0.85rem 0 0 1.25rem;
    white-space: pre;
    font-family: 'Prosto One', sans-serif;
    font-size: 1rem;
    line-height: 1.55;
    color: var(--accent);
    pointer-events: none;
    margin: 0;
  }

  /* Textarea styling */
  textarea {
    width: 100%;
    background: transparent;
    border: none;
    color: var(--text-primary);
    font-family: 'Prosto One', sans-serif;
    font-size: 1rem;
    line-height: 1.55;
    padding: 0.85rem 0.85rem 0.85rem 2.25rem;
    outline: none;
    resize: vertical;
    overflow: hidden;
  }

  textarea::placeholder {
    color: var(--text-tertiary);
  }

  textarea:disabled {
    opacity: 0.7;
    cursor: not-allowed;
  }

  @media (max-width: 640px) {
    .multi-input-wrapper {
      width: 100%;
      margin-bottom: 1.25rem;
    }
  }
</style>
