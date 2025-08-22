<script lang="ts">
  import { onMount } from 'svelte';
  import { DynamicIcon } from '$lib/components/icons';

  export let molecules = [
    { id: 'caffeine', name: 'Caffeine', formula: 'C₈H₁₀N₄O₂', mw: 194.19 },
    { id: 'aspirin', name: 'Aspirin', formula: 'C₉H₈O₄', mw: 180.16 },
    { id: 'glucose', name: 'Glucose', formula: 'C₆H₁₂O₆', mw: 180.16 },
  ];
  export const smiles: string = '';
  export const highlights: number[] = [];
  export const width: number = 300;
  export const height: number = 200;
  export const interactive: boolean = false;

  let selectedMolecule = molecules[0];
  let rotation = 0;
  let animationFrame: number;

  onMount(() => {
    const animate = () => {
      rotation += 0.5;
      animationFrame = requestAnimationFrame(animate);
    };
    animate();

    return () => {
      if (animationFrame) {
        cancelAnimationFrame(animationFrame);
      }
    };
  });
</script>

<div class="molecule-viewer">
  <div class="viewer-3d">
    <div class="molecule-container" style="transform: rotateY({rotation}deg)">
      <div class="molecule-structure">
        <div class="atom carbon" style="--x: 0; --y: 0; --z: 0"></div>
        <div class="atom oxygen" style="--x: 50; --y: 0; --z: 0"></div>
        <div class="atom nitrogen" style="--x: -50; --y: 0; --z: 0"></div>
        <div class="atom carbon" style="--x: 0; --y: 50; --z: 0"></div>
        <div class="atom carbon" style="--x: 0; --y: -50; --z: 0"></div>
        <div class="bond" style="--angle: 0deg"></div>
        <div class="bond" style="--angle: 72deg"></div>
        <div class="bond" style="--angle: 144deg"></div>
      </div>
    </div>
  </div>

  <div class="molecule-info">
    <h3>{selectedMolecule.name}</h3>
    <div class="formula">{selectedMolecule.formula}</div>
    <div class="mw">MW: {selectedMolecule.mw} Da</div>
  </div>

  <div class="molecule-selector">
    {#each molecules as mol}
      <button
        class="molecule-button"
        class:active={selectedMolecule.id === mol.id}
        on:click={() => (selectedMolecule = mol)}
      >
        {mol.name}
      </button>
    {/each}
  </div>
</div>

<style>
  .molecule-viewer {
    height: 100%;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    padding: 1rem;
  }

  .viewer-3d {
    width: 200px;
    height: 200px;
    position: relative;
    perspective: 800px;
    margin-bottom: 2rem;
  }

  .molecule-container {
    width: 100%;
    height: 100%;
    position: relative;
    transform-style: preserve-3d;
    transition: transform 0.1s linear;
  }

  .molecule-structure {
    position: absolute;
    width: 100%;
    height: 100%;
    transform-style: preserve-3d;
  }

  .atom {
    position: absolute;
    width: 30px;
    height: 30px;
    border-radius: 50%;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%)
      translate3d(calc(var(--x) * 1px), calc(var(--y) * 1px), calc(var(--z) * 1px));
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
  }

  .atom.carbon {
    background: #333;
  }

  .atom.oxygen {
    background: #ff4444;
  }

  .atom.nitrogen {
    background: #4444ff;
  }

  .bond {
    position: absolute;
    width: 60px;
    height: 4px;
    background: rgba(0, 0, 0, 0.3);
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%) rotate(var(--angle));
    transform-origin: center;
  }

  .molecule-info {
    text-align: center;
    margin-bottom: 1.5rem;
  }

  .molecule-info h3 {
    margin: 0 0 0.5rem;
    font-size: 1.25rem;
    color: var(--text-primary);
  }

  .formula {
    font-size: 1.1rem;
    color: var(--accent);
    margin-bottom: 0.25rem;
  }

  .mw {
    font-size: 0.9rem;
    color: var(--text-secondary);
  }

  .molecule-selector {
    display: flex;
    gap: 0.5rem;
  }

  .molecule-button {
    padding: 0.5rem 1rem;
    background: rgba(255, 255, 255, 0.6);
    border: 1px solid rgba(120, 121, 255, 0.3);
    border-radius: 20px;
    color: var(--text-primary);
    font-size: 0.85rem;
    cursor: pointer;
    transition: all 0.2s ease;
  }

  .molecule-button:hover {
    background: rgba(120, 121, 255, 0.1);
  }

  .molecule-button.active {
    background: var(--accent);
    color: white;
    border-color: var(--accent);
  }
</style>
