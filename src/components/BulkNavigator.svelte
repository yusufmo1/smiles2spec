<script>
  import { createEventDispatcher } from 'svelte';
  
  export let name = "";
  export let index = 0;   // 1-based
  export let total = 0;
  
  const dispatch = createEventDispatcher();
</script>

<div class="navigator glass-card">
  <div class="spacer"></div>
  
  <div class="center-group">
    <button class="arrow" on:click={() => dispatch('prev')} disabled={total<=1}>◀</button>
    <div class="name-container">
      <span class="name">{name || "—"}</span>
      {#if total}
        <span class="count">{index}/{total}</span>
      {/if}
    </div>
    <button class="arrow" on:click={() => dispatch('next')} disabled={total<=1}>▶</button>
  </div>
  
  <div class="spacer-right"></div>
</div>

<style>
.navigator {
  margin-top: -0.5rem;
  width: 100%;
  padding: 0.75rem 1.25rem;
  display: grid;
  grid-template-columns: 1fr auto 1fr;
  align-items: center;
  font-weight: 500;
  font-size: 0.95rem;
  border-radius: var(--enforce-pill);
}

.center-group {
  display: flex;
  align-items: center;
  gap: 0.75rem;
}

.spacer, .spacer-right {
  /* Empty spacer for grid layout */
}

.name-container {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  text-align: center;
}

.name {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  max-width: 250px;
}

.count { 
  color: var(--text-tertiary); 
  font-size: 0.8rem; 
}

.arrow {
  background: transparent;
  border: none;
  font-size: 1.1rem;
  color: var(--accent);
  cursor: pointer;
  padding: 0;
  width: 1.5rem;
  height: 1.5rem;
  display: flex;
  align-items: center;
  justify-content: center;
}

.arrow:disabled { 
  opacity: 0.25; 
  cursor: default; 
}
</style> 