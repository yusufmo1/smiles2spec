<script>
  import { createEventDispatcher } from 'svelte';
  
  // Track active page
  export let activePage = 'home';
  
  const dispatch = createEventDispatcher();
  
  // Handle navigation click
  function handleNavClick(page) {
    // Dispatch an event to notify the parent component
    dispatch('click', page);
    activePage = page;
  }
</script>

<nav class="navbar glass-card">
  <div class="navbar-container">
    <div class="logo">
      <span class="text-gradient prosto-one-regular">MS PREDICTOR</span>
    </div>
    <ul class="nav-links">
      <li>
        <button 
          class={activePage === 'home' ? 'active' : ''} 
          on:click={() => handleNavClick('home')}>
          Spectral Simulation
        </button>
      </li>
      <li>
        <button 
          class={activePage === 'how-it-works' ? 'active' : ''} 
          on:click={() => handleNavClick('how-it-works')}>
          How It Works
        </button>
      </li>
      <li>
        <button 
          class={activePage === 'chat' ? 'active' : ''} 
          on:click={() => handleNavClick('chat')}>
          Chat with Spectra
        </button>
      </li>
    </ul>
  </div>
</nav>

<style>
  .navbar {
    width: 100%;
    max-width: 90%;
    margin: 1.5rem auto 1rem;
    padding: 0.75rem 1.75rem;
    border-radius: var(--enforce-pill);
    display: flex;
    justify-content: center;
    transition: all var(--transition-smooth);
    backdrop-filter: blur(30px) saturate(160%);
    -webkit-backdrop-filter: blur(30px) saturate(160%);
    background: rgba(255,255,255,0.35);
    box-shadow:
      0 25px 55px rgba(0,0,0,0.08),
      0 3px 12px rgba(0,0,0,0.04);
    z-index: 10;
  }
  
  .navbar-container {
    width: 100%;
    display: flex;
    justify-content: space-between;
    align-items: center;
  }
  
  .logo {
    font-size: 1.25rem;
    font-weight: 600;
  }
  
  .nav-links {
    display: flex;
    list-style: none;
    margin: 0;
    padding: 0;
    gap: 1.5rem;
  }
  
  .nav-links button {
    background: transparent;
    border: none;
    color: var(--text-secondary);
    font-size: 1rem;
    font-weight: 500;
    cursor: pointer;
    padding: 0.5rem 1rem;
    border-radius: var(--enforce-pill);
    transition: all var(--transition-smooth);
  }
  
  .nav-links button:hover {
    color: var(--accent);
    background: rgba(120, 121, 255, 0.08);
  }
  
  .nav-links button.active {
    color: white;
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
    box-shadow: 0 4px 12px rgba(120, 121, 255, 0.2);
  }
  
  /* Responsive adjustments */
  @media (max-width: 640px) {
    .navbar {
      max-width: 95%;
      padding: 0.75rem 1.25rem;
      margin: 1rem auto 0.75rem;
    }
    
    .logo {
      font-size: 1.1rem;
    }
    
    .nav-links {
      gap: 0.75rem;
    }
    
    .nav-links button {
      font-size: 0.9rem;
      padding: 0.4rem 0.8rem;
    }
  }
  
  /* Extra-large screen adjustments */
  @media (min-width: 1920px) {
    .navbar {
      max-width: 85%;
    }
  }
</style> 