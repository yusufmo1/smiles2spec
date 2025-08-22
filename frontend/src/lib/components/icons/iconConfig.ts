export const IconSizes = {
  xs: 12,
  sm: 16,
  md: 20,
  lg: 24,
  xl: 32,
  '2xl': 48,
} as const;

export const IconColors = {
  primary: 'var(--text-primary)',
  secondary: 'var(--text-secondary)',
  tertiary: 'var(--text-tertiary)',
  accent: 'var(--accent)',
  success: '#34c759',
  error: '#ff453a',
  warning: '#ff9f0a',
} as const;

export type IconSize = keyof typeof IconSizes;
export type IconColor = keyof typeof IconColors;
