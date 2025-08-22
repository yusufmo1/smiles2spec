<script lang="ts">
  import { DynamicIcon } from '$lib/components/icons';
  import ProcessFlow from './components/ProcessFlow.svelte';

  const steps = [
    {
      id: 'validate',
      title: 'Validate',
      icon: 'ValidationIcon',
      color: '#34c759',
      description: 'Real-time SMILES syntax validation',
    },
    {
      id: 'parse',
      title: 'Parse',
      icon: 'MoleculeIcon',
      color: '#ff9f0a',
      description: 'Extract 8000+ molecular descriptors',
    },
    {
      id: 'predict',
      title: 'Predict',
      icon: 'CpuIcon',
      color: '#5856d6',
      description: 'AI predicts fragmentation patterns',
    },
    {
      id: 'generate',
      title: 'Generate',
      icon: 'SpectrumIcon',
      color: '#7879ff',
      description: 'Create interactive mass spectrum',
    },
  ];
</script>

<div class="workflow-panel">
  <ProcessFlow {steps} direction="horizontal" />
</div>

<style>
  .workflow-panel {
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
  }
</style>
