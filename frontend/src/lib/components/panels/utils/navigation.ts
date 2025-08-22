/**
 * Shared navigation utilities for panel components
 */

/**
 * Check if an element is interactive (button, input, etc.)
 */
export function isInteractiveElement(element: HTMLElement): boolean {
  // Only check for truly interactive elements that should block panel clicks
  const interactiveTags = ['BUTTON', 'TEXTAREA', 'INPUT', 'SELECT', 'A'];
  const interactiveClasses = ['send', 'export-btn', 'action-btn', 'cta-button'];

  // Direct tag check
  if (interactiveTags.includes(element.tagName)) {
    return true;
  }

  // Interactive attributes
  if (
    element.getAttribute('contenteditable') === 'true' ||
    element.getAttribute('role') === 'button' ||
    element.getAttribute('tabindex') === '0'
  ) {
    return true;
  }

  // Interactive classes
  if (interactiveClasses.some((cls) => element.classList.contains(cls))) {
    return true;
  }

  // Check if element is inside an actually interactive element
  if (
    element.closest(
      'button, textarea, input, select, a, [contenteditable="true"], .export-btn, .action-btn, .cta-button, .send'
    )
  ) {
    return true;
  }

  // Everything else is not interactive
  return false;
}

/**
 * Handle keyboard navigation for panels
 */
export function handlePanelKeyNavigation(
  event: KeyboardEvent,
  currentIndex: number,
  totalPanels: number,
  onNavigate: (index: number) => void,
  onClose?: () => void
): void {
  switch (event.key) {
    case 'ArrowLeft':
    case 'ArrowUp':
      event.preventDefault();
      const prevIndex = currentIndex > 0 ? currentIndex - 1 : totalPanels - 1;
      onNavigate(prevIndex);
      break;

    case 'ArrowRight':
    case 'ArrowDown':
      event.preventDefault();
      const nextIndex = currentIndex < totalPanels - 1 ? currentIndex + 1 : 0;
      onNavigate(nextIndex);
      break;

    case 'Home':
      event.preventDefault();
      onNavigate(0);
      break;

    case 'End':
      event.preventDefault();
      onNavigate(totalPanels - 1);
      break;

    case 'Escape':
      if (onClose) {
        event.preventDefault();
        onClose();
      }
      break;
  }
}

/**
 * Focus management utilities
 */
export class FocusManager {
  private previousFocus: HTMLElement | null = null;

  saveFocus(): void {
    this.previousFocus = document.activeElement as HTMLElement;
  }

  restoreFocus(): void {
    if (this.previousFocus && this.previousFocus.focus) {
      this.previousFocus.focus();
    }
  }

  trapFocus(container: HTMLElement): void {
    const focusableElements = container.querySelectorAll(
      'a[href], button:not([disabled]), textarea:not([disabled]), input:not([disabled]), select:not([disabled]), [tabindex]:not([tabindex="-1"])'
    );

    if (focusableElements.length === 0) return;

    const firstFocusable = focusableElements[0] as HTMLElement;
    const lastFocusable = focusableElements[focusableElements.length - 1] as HTMLElement;

    container.addEventListener('keydown', (e: KeyboardEvent) => {
      if (e.key !== 'Tab') return;

      if (e.shiftKey) {
        if (document.activeElement === firstFocusable) {
          e.preventDefault();
          lastFocusable.focus();
        }
      } else {
        if (document.activeElement === lastFocusable) {
          e.preventDefault();
          firstFocusable.focus();
        }
      }
    });

    // Focus first element
    setTimeout(() => firstFocusable.focus(), 100);
  }
}

/**
 * Debounce function for resize handlers
 */
export function debounce<T extends (...args: any[]) => void>(
  func: T,
  wait: number
): (...args: Parameters<T>) => void {
  let timeout: NodeJS.Timeout;

  return function executedFunction(...args: Parameters<T>) {
    const later = () => {
      clearTimeout(timeout);
      func(...args);
    };

    clearTimeout(timeout);
    timeout = setTimeout(later, wait);
  };
}
