import { writable, derived, get } from 'svelte/store';
import type { PageKey } from './index';
import { devLog } from '$lib/utils/performance';

// Panel definition interface
export interface Panel {
  id: string;
  title: string;
  component: string;
  props: Record<string, any>;
}

// Panel store state interface
export interface PanelStoreState {
  simulation: Panel[];
  info: Panel[];
  about: Panel[];
  chat: Panel[];
}

// Simplified panel definitions
const SIMULATION_PANELS: Panel[] = [
  { id: 'spectrum', title: 'MASS SPECTRUM', component: 'SpectrumPlot', props: {} },
  { id: 'fragmentation', title: 'TOP FRAGMENT IONS', component: 'IonFragmentation', props: {} },
  { id: 'structure', title: 'CHEMICAL STRUCTURE', component: 'StructurePanel', props: {} },
  { id: 'peaks', title: 'PEAK DATA TABLE', component: 'PeakTable', props: {} },
  { id: 'console', title: 'ANALYSIS CONSOLE', component: 'ConsoleOutput', props: {} },
  { id: 'export', title: 'EXPORT CENTER', component: 'ExportCenter', props: {} },
];

const INFO_PANELS: Panel[] = [
  {
    id: 'system-overview',
    title: 'SYSTEM OVERVIEW',
    component: 'SystemOverviewPanel',
    props: { compact: true },
  },
  {
    id: 'workflow',
    title: 'STEP-BY-STEP WORKFLOW',
    component: 'WorkflowPanel',
    props: { compact: true },
  },
  {
    id: 'ml-engine',
    title: 'MACHINE LEARNING ENGINE',
    component: 'MLEnginePanel',
    props: { compact: true },
  },
  {
    id: 'fragmentation-science',
    title: 'FRAGMENTATION SCIENCE',
    component: 'FragmentationSciencePanel',
    props: { compact: true },
  },
  {
    id: 'visualization-features',
    title: 'VISUALIZATION FEATURES',
    component: 'VisualizationFeaturesPanel',
    props: { compact: true },
  },
  {
    id: 'accuracy-validation',
    title: 'ACCURACY & VALIDATION',
    component: 'AccuracyValidationPanel',
    props: { compact: true },
  },
  {
    id: 'advanced-features',
    title: 'ADVANCED FEATURES',
    component: 'AdvancedFeaturesPanel',
    props: { compact: true },
  },
  {
    id: 'technical-requirements',
    title: 'TECHNICAL REQUIREMENTS',
    component: 'TechnicalRequirementsPanel',
    props: { compact: true },
  },
];

const ABOUT_PANELS: Panel[] = [
  {
    id: 'developer-hero',
    title: 'MEET THE DEVELOPER',
    component: 'DeveloperHeroPanel',
    props: { compact: true },
  },
  {
    id: 'developer-journey',
    title: 'MY JOURNEY',
    component: 'DeveloperJourneyPanel',
    props: { compact: true },
  },
  {
    id: 'project-metrics',
    title: 'PROJECT IMPACT',
    component: 'ProjectMetricsPanel',
    props: { compact: true },
  },
  { id: 'tech-stack', title: 'TECHNOLOGY', component: 'TechStackBanner', props: { compact: true } },
  {
    id: 'features-showcase',
    title: 'KEY FEATURES',
    component: 'FeaturesShowcasePanel',
    props: { compact: true },
  },
  {
    id: 'open-source',
    title: 'OPEN SOURCE',
    component: 'OpenSourceCommitmentPanel',
    props: { compact: true },
  },
  {
    id: 'contact-connect',
    title: 'CONNECT',
    component: 'ContactConnectRevampedPanel',
    props: { compact: true },
  },
  {
    id: 'acknowledgments',
    title: 'ACKNOWLEDGMENTS',
    component: 'AcknowledgmentsStreamlinedPanel',
    props: { compact: true },
  },
];

const CHAT_PANELS: Panel[] = [
  {
    id: 'chat-conversation',
    title: 'CHAT WITH SPECTRUM',
    component: 'ChatWithSpectrum',
    props: {},
  },
  { id: 'chat-spectrum', title: 'PREDICTED SPECTRUM', component: 'SpectrumPlot', props: {} },
  {
    id: 'compound-stats',
    title: 'WHAT SPECTRUM KNOWS',
    component: 'CompoundStatsPanel',
    props: {},
  },
];

function createPanelStore() {
  const initialState: PanelStoreState = {
    simulation: SIMULATION_PANELS.map((p) => ({ ...p, props: { ...p.props } })),
    info: INFO_PANELS.map((p) => ({ ...p, props: { ...p.props } })),
    about: ABOUT_PANELS.map((p) => ({ ...p, props: { ...p.props } })),
    chat: CHAT_PANELS.map((p) => ({ ...p, props: { ...p.props } })),
  };

  const { subscribe, set, update } = writable<PanelStoreState>(initialState);

  return {
    subscribe,
    updatePanelProps: (panelId: string, newProps: Record<string, any>) => {
      devLog.log(`Updating panel ${panelId} with props:`, newProps);

      update((state) => {
        const newState: PanelStoreState = {
          ...state,
          simulation: state.simulation.map((panel) => {
            if (panel.id === panelId) {
              const updatedPanel = {
                ...panel,
                props: { ...panel.props, ...newProps },
              };
              devLog.log(`Panel ${panelId} updated:`, updatedPanel.props);
              return updatedPanel;
            }
            return panel;
          }),
          // Also update other categories in case the panel is there
          info: state.info.map((panel) => {
            if (panel.id === panelId) {
              const updatedPanel = {
                ...panel,
                props: { ...panel.props, ...newProps },
              };
              devLog.log(`Panel ${panelId} updated in info:`, updatedPanel.props);
              return updatedPanel;
            }
            return panel;
          }),
          about: state.about.map((panel) => {
            if (panel.id === panelId) {
              const updatedPanel = {
                ...panel,
                props: { ...panel.props, ...newProps },
              };
              devLog.log(`Panel ${panelId} updated in about:`, updatedPanel.props);
              return updatedPanel;
            }
            return panel;
          }),
          chat: state.chat.map((panel) => {
            if (panel.id === panelId) {
              const updatedPanel = {
                ...panel,
                props: { ...panel.props, ...newProps },
              };
              devLog.log(`Panel ${panelId} updated in chat:`, updatedPanel.props);
              return updatedPanel;
            }
            return panel;
          }),
        };

        devLog.log(`Store state updated for ${panelId}`);
        return newState;
      });
    },
    getPanelsForPage: (page: PageKey): Panel[] => {
      const currentState = get({ subscribe });
      switch (page) {
        case 'home':
          return []; // Landing page has no panels
        case 'spectral-simulation':
          return currentState.simulation;
        case 'how-it-works':
          return currentState.info;
        case 'about':
          return currentState.about;
        case 'chat-with-spectrum':
          return currentState.chat;
        default:
          return currentState.simulation;
      }
    },
  };
}

export const panelStore = createPanelStore();

// Enhanced derived stores with debugging
export const simulationPanels = derived([panelStore], ([$store]) => {
  devLog.log('simulationPanels derived store updated');
  return $store.simulation;
});

export const infoPanels = derived([panelStore], ([$store]) => {
  devLog.log('infoPanels derived store updated');
  return $store.info;
});

export const aboutPanels = derived([panelStore], ([$store]) => {
  devLog.log('aboutPanels derived store updated');
  return $store.about;
});

export const chatPanels = derived([panelStore], ([$store]) => {
  devLog.log('chatPanels derived store updated');
  return $store.chat;
});
