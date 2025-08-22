<script>
  import IconShowcase from '$lib/components/debug/IconShowcase.svelte';
  import { DynamicIcon } from '$lib/components/icons';
</script>

<svelte:head>
  <title>Icon System - Mass Spectrum Predictor</title>
</svelte:head>

<div class="icon-docs">
  <header>
    <h1>Icon System Documentation</h1>
    <p>A comprehensive guide to using the Mass Spectrum Predictor icon system</p>
  </header>

  <section class="usage-guide">
    <h2>Usage Guide</h2>

    <div class="code-sample">
      <h3>Direct Component Import</h3>
      <pre>import {DynamicIcon} from '$lib/components/icons';

&lt;DynamicIcon name="CheckIcon" size=&#123;24&#125; color="var(--accent)" /&gt;</pre>
    </div>

    <div class="code-sample">
      <h3>Using the DynamicIcon Component</h3>
      <pre>import {DynamicIcon} from '$lib/components/icons';

&lt;DynamicIcon name="CheckIcon" size=&#123;24&#125; color="var(--accent)" /&gt;</pre>
    </div>

    <div class="code-sample">
      <h3>Dynamic Icon with Props</h3>
      <pre>import {DynamicIcon} from '$lib/components/icons';

&lt;DynamicIcon 
  name="SpectrumIcon" 
  size=&#123;32&#125; 
  color="var(--accent)" 
  aria-label="Mass spectrum visualization" 
/&gt;</pre>
    </div>
  </section>

  <section class="accessibility">
    <h2>Accessibility</h2>
    <p>Icons should include ARIA attributes when they convey meaning:</p>

    <pre>&lt;DynamicIcon 
  name="SpectrumIcon"
  size=&#123;24&#125; 
  color="var(--text-primary)" 
  aria-label="Chart data" 
  role="img"
/&gt;</pre>
  </section>

  <section>
    <h2>Icon Library</h2>
    <IconShowcase />
  </section>
</div>

<style>
  .icon-docs {
    max-width: 1200px;
    margin: 0 auto;
    padding: 2rem;
  }

  header {
    margin-bottom: 3rem;
    text-align: center;
  }

  h1 {
    margin: 0 0 0.5rem;
    font-size: 2.5rem;
    color: var(--text-primary);
  }

  header p {
    font-size: 1.2rem;
    color: var(--text-secondary);
    margin: 0;
    opacity: 0.8;
  }

  section {
    margin-bottom: 3rem;
  }

  h2 {
    margin: 0 0 1.5rem;
    font-size: 1.75rem;
    color: var(--text-primary);
    border-bottom: 1px solid var(--surface-stroke);
    padding-bottom: 0.5rem;
  }

  h3 {
    margin: 1.5rem 0 0.75rem;
    font-size: 1.25rem;
    color: var(--text-primary);
  }

  .code-sample {
    margin-bottom: 2rem;
  }

  pre {
    background: rgba(0, 0, 0, 0.03);
    padding: 1.5rem;
    border-radius: 8px;
    overflow-x: auto;
    font-family: monospace;
    font-size: 0.9rem;
    line-height: 1.5;
    margin: 0.75rem 0;
  }

  .accessibility {
    background: rgba(120, 121, 255, 0.05);
    padding: 1.5rem;
    border-radius: 8px;
    border: 1px solid var(--accent-soft);
    margin-bottom: 3rem;
  }

  .accessibility h2 {
    margin-top: 0;
    border-bottom-color: rgba(120, 121, 255, 0.2);
  }
</style>
