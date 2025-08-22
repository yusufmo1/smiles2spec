// Icon loader with explicit imports for production compatibility
import type { ComponentType } from 'svelte';

export type IconName =
  | 'GridIcon'
  | 'CarouselIcon'
  | 'CloseIcon'
  | 'ArrowLeftIcon'
  | 'ArrowRightIcon'
  | 'SearchIcon'
  | 'UploadIcon'
  | 'DownloadIcon'
  | 'GenerateIcon'
  | 'SendIcon'
  | 'CopyIcon'
  | 'ExportIcon'
  | 'MoleculeIcon'
  | 'SpectrumIcon'
  | 'CpuIcon'
  | 'DatabaseIcon'
  | 'FragmentationIcon'
  | 'DataTableIcon'
  | 'ChatWithSpectrumIcon'
  | 'CheckIcon'
  | 'ErrorIcon'
  | 'InfoIcon'
  | 'LoadingIcon'
  | 'WandIcon'
  | 'FlaskIcon'
  | 'SvelteKitIcon'
  | 'RDKitIcon'
  | 'ScikitLearnIcon'
  | 'DockerIcon'
  | 'TypeScriptIcon'
  | 'CSSIcon'
  | 'PlotlyJSIcon'
  | 'PythonIcon'
  | 'NPMIcon'
  | 'MatplotlibIcon'
  | 'PyTorchIcon'
  | 'TensorFlowIcon'
  | 'UserIcon'
  | 'GraduationCapIcon'
  | 'MicroscopeIcon'
  | 'LightbulbIcon'
  | 'SettingsIcon'
  | 'BookIcon'
  | 'HeartIcon'
  | 'StarIcon'
  | 'HomeIcon'
  | 'PhoneIcon'
  | 'MailIcon'
  | 'TwitterIcon'
  | 'LinkedinIcon'
  | 'GithubIcon'
  | 'AtomIcon'
  | 'BondIcon'
  | 'TestTubeIcon'
  | 'ValidationIcon'
  | 'RefreshIcon'
  | 'FlashIcon'
  | 'BeakerIcon'
  | 'TargetIcon'
  | 'ScaleIcon'
  | 'ChartIcon'
  | 'MenuIcon';

// Cache for loaded icons
const iconCache = new Map<IconName, ComponentType>();

// Icon map with explicit imports for production compatibility
export const ICON_MAP: Record<IconName, () => Promise<ComponentType>> = {
  GridIcon: () => import('./GridIcon.svelte').then((m) => m.default),
  CarouselIcon: () => import('./CarouselIcon.svelte').then((m) => m.default),
  CloseIcon: () => import('./CloseIcon.svelte').then((m) => m.default),
  ArrowLeftIcon: () => import('./ArrowLeftIcon.svelte').then((m) => m.default),
  ArrowRightIcon: () => import('./ArrowRightIcon.svelte').then((m) => m.default),
  SearchIcon: () => import('./SearchIcon.svelte').then((m) => m.default),
  UploadIcon: () => import('./UploadIcon.svelte').then((m) => m.default),
  DownloadIcon: () => import('./DownloadIcon.svelte').then((m) => m.default),
  GenerateIcon: () => import('./GenerateIcon.svelte').then((m) => m.default),
  SendIcon: () => import('./SendIcon.svelte').then((m) => m.default),
  CopyIcon: () => import('./CopyIcon.svelte').then((m) => m.default),
  ExportIcon: () => import('./ExportIcon.svelte').then((m) => m.default),
  MoleculeIcon: () => import('./MoleculeIcon.svelte').then((m) => m.default),
  SpectrumIcon: () => import('./SpectrumIcon.svelte').then((m) => m.default),
  CpuIcon: () => import('./CpuIcon.svelte').then((m) => m.default),
  DatabaseIcon: () => import('./DatabaseIcon.svelte').then((m) => m.default),
  FragmentationIcon: () => import('./FragmentationIcon.svelte').then((m) => m.default),
  DataTableIcon: () => import('./DataTableIcon.svelte').then((m) => m.default),
  ChatWithSpectrumIcon: () => import('./ChatWithSpectrumIcon.svelte').then((m) => m.default),
  CheckIcon: () => import('./CheckIcon.svelte').then((m) => m.default),
  ErrorIcon: () => import('./ErrorIcon.svelte').then((m) => m.default),
  InfoIcon: () => import('./InfoIcon.svelte').then((m) => m.default),
  LoadingIcon: () => import('./LoadingIcon.svelte').then((m) => m.default),
  WandIcon: () => import('./WandIcon.svelte').then((m) => m.default),
  FlaskIcon: () => import('./FlaskIcon.svelte').then((m) => m.default),
  SvelteKitIcon: () => import('./SvelteKitIcon.svelte').then((m) => m.default),
  RDKitIcon: () => import('./RDKitIcon.svelte').then((m) => m.default),
  ScikitLearnIcon: () => import('./ScikitLearnIcon.svelte').then((m) => m.default),
  DockerIcon: () => import('./DockerIcon.svelte').then((m) => m.default),
  TypeScriptIcon: () => import('./TypeScriptIcon.svelte').then((m) => m.default),
  CSSIcon: () => import('./CSSIcon.svelte').then((m) => m.default),
  PlotlyJSIcon: () => import('./PlotlyJSIcon.svelte').then((m) => m.default),
  PythonIcon: () => import('./PythonIcon.svelte').then((m) => m.default),
  NPMIcon: () => import('./NPMIcon.svelte').then((m) => m.default),
  MatplotlibIcon: () => import('./MatplotlibIcon.svelte').then((m) => m.default),
  PyTorchIcon: () => import('./PyTorchIcon.svelte').then((m) => m.default),
  TensorFlowIcon: () => import('./TensorFlowIcon.svelte').then((m) => m.default),
  UserIcon: () => import('./UserIcon.svelte').then((m) => m.default),
  GraduationCapIcon: () => import('./GraduationCapIcon.svelte').then((m) => m.default),
  MicroscopeIcon: () => import('./MicroscopeIcon.svelte').then((m) => m.default),
  LightbulbIcon: () => import('./LightbulbIcon.svelte').then((m) => m.default),
  SettingsIcon: () => import('./SettingsIcon.svelte').then((m) => m.default),
  BookIcon: () => import('./BookIcon.svelte').then((m) => m.default),
  HeartIcon: () => import('./HeartIcon.svelte').then((m) => m.default),
  StarIcon: () => import('./StarIcon.svelte').then((m) => m.default),
  HomeIcon: () => import('./HomeIcon.svelte').then((m) => m.default),
  PhoneIcon: () => import('./PhoneIcon.svelte').then((m) => m.default),
  MailIcon: () => import('./MailIcon.svelte').then((m) => m.default),
  TwitterIcon: () => import('./TwitterIcon.svelte').then((m) => m.default),
  LinkedinIcon: () => import('./LinkedinIcon.svelte').then((m) => m.default),
  GithubIcon: () => import('./GithubIcon.svelte').then((m) => m.default),
  AtomIcon: () => import('./AtomIcon.svelte').then((m) => m.default),
  BondIcon: () => import('./BondIcon.svelte').then((m) => m.default),
  TestTubeIcon: () => import('./TestTubeIcon.svelte').then((m) => m.default),
  ValidationIcon: () => import('./ValidationIcon.svelte').then((m) => m.default),
  RefreshIcon: () => import('./RefreshIcon.svelte').then((m) => m.default),
  FlashIcon: () => import('./FlashIcon.svelte').then((m) => m.default),
  BeakerIcon: () => import('./BeakerIcon.svelte').then((m) => m.default),
  TargetIcon: () => import('./TargetIcon.svelte').then((m) => m.default),
  ScaleIcon: () => import('./ScaleIcon.svelte').then((m) => m.default),
  ChartIcon: () => import('./ChartIcon.svelte').then((m) => m.default),
  MenuIcon: () => import('./MenuIcon.svelte').then((m) => m.default),
};

// Helper function for using the icon map with caching
export async function loadIcon(iconName: IconName): Promise<ComponentType> {
  if (iconCache.has(iconName)) {
    return iconCache.get(iconName)!;
  }

  const iconLoader = ICON_MAP[iconName];
  if (!iconLoader) {
    throw new Error(`Icon "${iconName}" not found in ICON_MAP`);
  }

  const IconComponent = await iconLoader();
  iconCache.set(iconName, IconComponent);
  return IconComponent;
}

export { default as Icon } from './Icon.svelte';
export { default as DynamicIcon } from './DynamicIcon.svelte';
export { IconSizes, IconColors, type IconSize, type IconColor } from './iconConfig';

export interface IconProps {
  size?: number | string;
  color?: string;
  className?: string;
}

// ===== ICONS FIX COMPLETE =====
// Icons now use explicit imports in ICON_MAP for production compatibility
// Use: <DynamicIcon name="IconName" /> for lazy loading with production support
