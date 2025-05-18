<script>
  export let className = "";
</script>

<header>
  <div class="title-container">
    <h1 class="display text-gradient {className}">SPECTRAL SIMULATION</h1>
    <p>Enter SMILES string to predict high-resolution LC-ESI-qTOF mass spectra</p>
  </div>
</header>

<style>
  header {
    margin-bottom: 2rem;
    text-align: center;
    padding: 0;
    position: relative;
    overflow: hidden;
    height: auto;
  }
  
  .title-container {
    position: relative;
    z-index: 2;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
  }
  
  .display {
    font: 600 clamp(2.75rem, 5vw, 4.25rem)/1 'SF Pro Display', system-ui;
    margin: 0 0 1rem;
    letter-spacing: 0.05em;
  }
  
  p {
    margin: 0;
    color: var(--text-secondary);
    font-size: 1.1rem;
    font-weight: 400;
  }
  
  .text-gradient {
    background: linear-gradient(135deg, #7879ff 0%, #6d56f7 60%, #7879ff 100%) 0 0/200%;
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    text-fill-color: transparent;
    animation: title-shift 15s linear infinite;
  }
  
  @keyframes title-shift {
    to { background-position: 200% 0; }
  }
  
  :global(.prosto-one-regular) {
    font-family: "Prosto One", sans-serif !important;
    font-weight: 400;
    font-style: normal;
  }
</style> 