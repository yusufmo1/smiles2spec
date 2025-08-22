import type { IconName } from '$lib/components/icons';

export type IconType =
  | 'error'
  | 'success'
  | 'info'
  | 'loading'
  | 'cpu'
  | 'chart'
  | 'generate'
  | 'flask'
  | 'molecule'
  | 'spectrum'
  | 'database'
  | 'fragmentation'
  | 'datatable'
  | 'chat'
  | 'microscope'
  | 'refresh'
  | 'settings'
  | 'atom'
  | 'bond'
  | 'beaker';

/**
 * Maps semantic icon types to their corresponding icon names for DynamicIcon
 * @param type Semantic icon type
 * @returns Icon name string for DynamicIcon component
 */
export const getIconByType = (type: IconType): IconName => {
  const iconMap: Record<IconType, IconName> = {
    error: 'ErrorIcon',
    success: 'CheckIcon',
    info: 'InfoIcon',
    loading: 'LoadingIcon',
    cpu: 'CpuIcon',
    chart: 'SpectrumIcon',
    generate: 'WandIcon',
    flask: 'FlaskIcon',
    molecule: 'MoleculeIcon',
    spectrum: 'SpectrumIcon',
    database: 'DatabaseIcon',
    fragmentation: 'FragmentationIcon',
    datatable: 'DataTableIcon',
    chat: 'ChatWithSpectrumIcon',
    microscope: 'MicroscopeIcon',
    refresh: 'RefreshIcon',
    settings: 'SettingsIcon',
    atom: 'AtomIcon',
    bond: 'BondIcon',
    beaker: 'BeakerIcon',
  };

  return iconMap[type] || 'InfoIcon';
};

/**
 * Returns the ARIA label for a given icon type
 * @param type Semantic icon type
 * @returns Appropriate ARIA label
 */
export const getIconAriaLabel = (type: IconType): string => {
  const ariaLabels: Record<IconType, string> = {
    error: 'Error',
    success: 'Success',
    info: 'Information',
    loading: 'Loading',
    cpu: 'CPU',
    chart: 'Chart',
    generate: 'Generate',
    flask: 'Laboratory flask',
    molecule: 'Molecule',
    spectrum: 'Spectrum',
    database: 'Database',
    fragmentation: 'Fragmentation',
    datatable: 'Data table',
    chat: 'Chat with spectrum',
    microscope: 'Microscope',
    refresh: 'Refresh',
    settings: 'Settings',
    atom: 'Atom',
    bond: 'Bond',
    beaker: 'Beaker',
  };

  return ariaLabels[type] || '';
};
