import { writable, derived } from 'svelte/store';
import { carouselIndex } from '../stores.js';

export const panelDefinitions = writable([
  { id: 'spectrum', title: 'MASS SPECTRUM', component: 'SpectrumPlot', props: {} },
  { id: 'fragments', title: 'TOP FRAGMENT IONS', component: 'IonFragmentation', props: {} },
  { id: 'structure', title: 'CHEMICAL STRUCTURE', component: 'StructurePanel', props: {} },
  { id: 'peaks', title: 'PEAK DATA TABLE', component: 'PeakTable', props: {} },
  { id: 'console', title: 'ANALYSIS CONSOLE', component: 'ConsoleOutput', props: {} },
  { id: 'chat', title: 'CHAT WITH SPECTRUM', component: 'ChatWithSpectrum', props: {} }
]);

export const updatePanelProps = (id, props) => {
  panelDefinitions.update(panels =>
    panels.map(panel =>
      panel.id === id ? { ...panel, props: { ...panel.props, ...props } } : panel
    )
  );
};

export const activePanelData = derived(
  [panelDefinitions, carouselIndex],
  ([$panels, $index]) => $panels[$index] || null
);
