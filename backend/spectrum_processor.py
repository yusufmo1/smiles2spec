"""
Spectral data processing and binning functionality.
"""

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from .utils import handle_error, logger, ensure_numpy_array

class SpectralProcessor:
    def __init__(self, config):
        self.config = config
    
    def process_peaks(self, peaks_data):
        """Process a list of peaks (mz, intensity) into a standardized format."""
        try:
            # Initialize required result variables
            max_peaks = self.config.get('max_peaks', 499)
            peak_array = np.zeros((max_peaks, 2), dtype=np.float32)
            attention_mask = np.zeros(max_peaks, dtype=np.float32)
            
            # Convert to DataFrame if list of dicts
            if isinstance(peaks_data, list) and all(isinstance(p, dict) for p in peaks_data):
                df = pd.DataFrame(peaks_data)
            elif isinstance(peaks_data, list) and all(isinstance(p, (list, tuple)) for p in peaks_data):
                df = pd.DataFrame(peaks_data, columns=["mz", "intensity"])
            elif isinstance(peaks_data, pd.DataFrame):
                df = peaks_data
            else:
                raise ValueError("Peaks data must be a DataFrame or list of dicts/tuples with mz and intensity")
            
            # Ensure we have the right columns
            if 'mz' not in df.columns or 'intensity' not in df.columns:
                raise ValueError("Peaks data must have 'mz' and 'intensity' columns")
            
            # Apply binning if configured
            if self.config.get('bin_size') and not df.empty:
                df['binned_mz'] = (df['mz'] / self.config['bin_size']).astype(int) * self.config['bin_size']
                df = df.groupby('binned_mz').agg({'intensity': 'max'}).reset_index().rename(columns={'binned_mz': 'mz'})
            
            # Process peaks if data is not empty
            if not df.empty:
                sort_by = self.config.get('sort_peaks_by', 'intensity')
                sorted_peaks = df.sort_values(sort_by, ascending=False)
                max_intensity = sorted_peaks['intensity'].max() if not sorted_peaks.empty else 0
                
                if max_intensity > 0:
                    sorted_peaks['intensity'] /= max_intensity
                
                peak_sequence = list(zip(sorted_peaks['mz'], sorted_peaks['intensity']))
                original_peak_count = len(peak_sequence)
                
                if len(peak_sequence) > 0:
                    # Process peak sequences
                    if len(peak_sequence) > max_peaks:
                        attention_mask = np.ones(max_peaks)
                        peak_sequence = peak_sequence[:max_peaks]
                    else:
                        attention_mask = np.concatenate([
                            np.ones(len(peak_sequence)), 
                            np.zeros(max_peaks - len(peak_sequence))
                        ])
                        peak_sequence += [(0.0, 0.0)] * (max_peaks - len(peak_sequence))
                    
                    peak_array = np.array(peak_sequence, dtype=np.float32)
            
            # Get the true maximum m/z value
            true_max_mz = df['mz'].max() if not df.empty else 0
            
            return {
                'peaks': peak_array,
                'attention_mask': attention_mask,
                'original_peak_count': len(df) if not df.empty else 0,
                'max_mz': true_max_mz,
                'true_max_mz': true_max_mz,
                'bin_size': self.config.get('bin_size', 1.0)
            }
        except Exception as e:
            error_result = handle_error(e, "processing peaks")
            logger.error(f"Error processing peaks: {error_result['error']}")
            raise

    def create_binned_spectrum(self, peaks_data):
        """Create a binned spectrum from peaks data."""
        try:
            # Process the peaks to get standardized format
            processed = self.process_peaks(peaks_data)
            peaks = processed['peaks']
            mask = processed['attention_mask']
            real_peaks = peaks[mask.astype(bool)]
            
            # Create binned spectrum up to configured max_mz
            max_mz = self.config.get('max_mz', 499)
            bin_size = self.config.get('bin_size', 1.0)
            num_bins = int(max_mz / bin_size) + 1
            binned_spectrum = np.zeros(num_bins)
            
            for mz, intensity in real_peaks:
                if mz < max_mz:
                    bin_idx = min(int(mz / bin_size), num_bins - 1)
                    binned_spectrum[bin_idx] = max(binned_spectrum[bin_idx], intensity)
            
            return {
                'binned_spectrum': binned_spectrum,
                'mz_values': np.arange(0, num_bins * bin_size, bin_size),
                'original_peaks': processed
            }
        except Exception as e:
            error_result = handle_error(e, "creating binned spectrum")
            logger.error(f"Error creating binned spectrum: {error_result['error']}")
            raise