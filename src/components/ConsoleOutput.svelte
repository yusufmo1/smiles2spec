<script>
  import { onMount, afterUpdate } from 'svelte';
  
  export let output = "Initializing...";
  let consoleElement;
  
  // Auto-scroll to bottom when content changes
  afterUpdate(() => {
    if (consoleElement) {
      consoleElement.scrollTop = consoleElement.scrollHeight;
    }
  });
</script>

<div class="console-container">
  <div class="console-content" bind:this={consoleElement}>
    <pre>{output}</pre>
  </div>
</div>

<style>
  .console-container {
    width: 100%;
    height: 100%;
    position: relative;
    min-height: 248px;
    border-radius: var(--enforce-pill);
    overflow: hidden;
  }
  
  .console-content {
    position: relative;
    height: 100%;
    max-height: 220px;
    padding: 1.25rem;
    background: rgba(255, 255, 255, 0.8);
    border-radius: var(--enforce-pill);
    overflow-y: auto;
    scrollbar-width: thin;
    font-family: 'SF Mono', 'Menlo', monospace;
    font-size: 0.8rem;
    line-height: 1.5;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
    border: 1px solid rgba(0, 0, 0, 0.05);
  }
  
  pre {
    margin: 0;
    white-space: pre-wrap;
    word-break: break-word;
    color: rgba(0, 0, 0, 0.85);
    overflow: visible;
  }
  
  /* Custom scrollbar for webkit browsers */
  .console-content::-webkit-scrollbar {
    width: 4px;
  }
  
  .console-content::-webkit-scrollbar-thumb {
    background-color: rgba(0, 0, 0, 0.15);
    border-radius: 2px;
  }
  
  .console-content::-webkit-scrollbar-track {
    background-color: rgba(0, 0, 0, 0.05);
  }
</style> 