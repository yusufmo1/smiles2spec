"""Spectrum processing functionality."""
import numpy as np
import pandas as pd
from typing import List, Dict, Any, Tuple
from scipy.signal import savgol_filter

from ...config import settings
from ...core.models.spectrum import Spectrum, Peak
from ...utils import logger

class SpectralProcessor:
    """
    Processor for mass spectral data operations.
    
    Handles conversion between different spectral representations,
    peak processing, binning, and normalization. Designed to work
    with both predicted and experimental mass spectra.
    
    Key operations:
    - Convert ML predictions to spectrum objects
    - Process peak lists with binning and filtering
    - Create standardized representations for visualization
    - Handle empty or invalid spectral data gracefully
    
    Attributes:
        config: Spectral processing configuration
        logger: Logger instance for error tracking
    """
    
    def __init__(self):
        self.config = settings.spectral
        self.logger = logger
    
    def create_spectrum_from_prediction(self, predicted_intensities: np.ndarray) -> Spectrum:
        """
        Create a Spectrum object from ML model predictions.
        
        Converts raw intensity predictions into a structured spectrum with
        appropriate m/z values based on the configured bin size.
        
        Args:
            predicted_intensities: Array of predicted intensities from ML model
            
        Returns:
            Spectrum object containing:
                - x: m/z values array
                - y: intensity values array
                - peaks: List of Peak objects for non-zero intensities
                
        Note:
            m/z values are generated based on bin_size configuration
            (default: 1 Da bins)
        """
        try:
            # Create m/z values based on bin size
            mz_values = np.arange(
                0, 
                len(predicted_intensities) * self.config.bin_size, 
                self.config.bin_size
            )
            
            # Ensure arrays have same length
            if len(mz_values) != len(predicted_intensities):
                min_len = min(len(mz_values), len(predicted_intensities))
                mz_values = mz_values[:min_len]
                predicted_intensities = predicted_intensities[:min_len]
            
            # Create peaks from non-zero intensities
            peaks = self._create_peaks_from_arrays(mz_values, predicted_intensities)
            
            return Spectrum(
                x=mz_values,
                y=predicted_intensities,
                peaks=peaks
            )
            
        except Exception as e:
            self.logger.error(f"Error creating spectrum from prediction: {str(e)}")
            raise
    
    def create_empty_spectrum(self) -> Spectrum:
        """
        Create an empty spectrum with zero intensities.
        
        Used as a fallback when no valid spectrum data is available.
        The spectrum covers the full m/z range defined in configuration.
        
        Returns:
            Spectrum object with:
                - x: Full m/z range array
                - y: Zero intensity array
                - peaks: Empty peak list
        """
        try:
            num_bins = int(self.config.max_mz / self.config.bin_size) + 1
            mz_values = np.arange(0, num_bins * self.config.bin_size, self.config.bin_size)
            intensities = np.zeros(num_bins)
            
            return Spectrum(
                x=mz_values,
                y=intensities,
                peaks=[]
            )
            
        except Exception as e:
            self.logger.error(f"Error creating empty spectrum: {str(e)}")
            raise
    
    def create_spectrum_from_peaks(self, peaks_data: List[Dict[str, float]]) -> Spectrum:
        """
        Create a binned spectrum from a list of peaks.
        
        Processes raw peak data into a binned representation suitable for
        visualization and comparison. Handles overlapping peaks by taking
        the maximum intensity per bin.
        
        Args:
            peaks_data: List of dictionaries with 'mz' and 'intensity' keys
            
        Returns:
            Spectrum object with binned representation
            
        Process:
            1. Bin peaks according to bin_size
            2. Aggregate overlapping peaks (max intensity)
            3. Create continuous spectrum representation
            4. Extract significant peaks for peak list
        """
        try:
            if not peaks_data:
                return self.create_empty_spectrum()
            
            # Convert to DataFrame for processing
            df = pd.DataFrame(peaks_data)
            
            # Ensure we have the right columns
            if 'mz' not in df.columns or 'intensity' not in df.columns:
                raise ValueError("Peaks data must have 'mz' and 'intensity' columns")
            
            # Apply binning if configured
            if self.config.bin_size and not df.empty:
                df['binned_mz'] = (df['mz'] / self.config.bin_size).astype(int) * self.config.bin_size
                df = df.groupby('binned_mz').agg({'intensity': 'max'}).reset_index()
                df = df.rename(columns={'binned_mz': 'mz'})
            
            # Create binned spectrum
            num_bins = int(self.config.max_mz / self.config.bin_size) + 1
            mz_values = np.arange(0, num_bins * self.config.bin_size, self.config.bin_size)
            binned_intensities = np.zeros(num_bins)
            
            # Fill bins with peak intensities
            for _, row in df.iterrows():
                mz, intensity = row['mz'], row['intensity']
                if mz < self.config.max_mz:
                    bin_idx = min(int(mz / self.config.bin_size), num_bins - 1)
                    binned_intensities[bin_idx] = max(binned_intensities[bin_idx], intensity)
            
            # Create peaks from processed data
            peaks = self._create_peaks_from_dataframe(df)
            
            return Spectrum(
                x=mz_values,
                y=binned_intensities,
                peaks=peaks
            )
            
        except Exception as e:
            self.logger.error(f"Error creating spectrum from peaks: {str(e)}")
            raise
    
    def process_peaks(self, peaks_data: List[Dict[str, float]]) -> Dict[str, Any]:
        """
        Process raw peaks into standardized format for ML models.
        
        Performs normalization, sorting, and padding to create consistent
        input format for models that work with peak lists rather than
        continuous spectra.
        
        Args:
            peaks_data: List of peak dictionaries
            
        Returns:
            Dictionary containing:
                - peaks: Padded array of [m/z, intensity] pairs
                - attention_mask: Binary mask for valid peaks
                - original_peak_count: Number of input peaks
                - processed_peak_count: Number after filtering
                - max_mz: Maximum m/z value in spectrum
                
        Processing Steps:
            1. Sort by intensity or m/z (configurable)
            2. Normalize intensities to [0, 1]
            3. Limit to max_peaks threshold
            4. Pad with zeros for consistent dimensions
        """
        try:
            if not peaks_data:
                return self._create_empty_peaks_result()
            
            # Convert to DataFrame
            df = pd.DataFrame(peaks_data)
            
            # Validate columns
            if 'mz' not in df.columns or 'intensity' not in df.columns:
                raise ValueError("Peaks data must have 'mz' and 'intensity' columns")
            
            # Sort by configured criterion
            sort_by = self.config.sort_peaks_by
            if sort_by in df.columns:
                df = df.sort_values(sort_by, ascending=False)
            
            # Normalize intensities
            if not df.empty:
                max_intensity = df['intensity'].max()
                if max_intensity > 0:
                    df['intensity'] = df['intensity'] / max_intensity
            
            # Limit to max peaks
            max_peaks = self.config.max_peaks
            if len(df) > max_peaks:
                df = df.head(max_peaks)
            
            # Create attention mask
            attention_mask = np.concatenate([
                np.ones(len(df)),
                np.zeros(max_peaks - len(df))
            ]) if len(df) < max_peaks else np.ones(max_peaks)
            
            # Create peak array
            peak_array = np.zeros((max_peaks, 2), dtype=np.float32)
            for i, (_, row) in enumerate(df.iterrows()):
                if i < max_peaks:
                    peak_array[i] = [row['mz'], row['intensity']]
            
            return {
                'peaks': peak_array,
                'attention_mask': attention_mask,
                'original_peak_count': len(peaks_data),
                'processed_peak_count': len(df),
                'max_mz': df['mz'].max() if not df.empty else 0
            }
            
        except Exception as e:
            self.logger.error(f"Error processing peaks: {str(e)}")
            raise
    
    def _create_peaks_from_arrays(self, mz_values: np.ndarray, intensities: np.ndarray, 
                                 threshold: float = 0.01) -> List[Peak]:
        """
        Create Peak objects from parallel arrays of m/z and intensity.
        
        Args:
            mz_values: Array of m/z values
            intensities: Array of intensity values
            threshold: Minimum intensity threshold (relative to max)
            
        Returns:
            List of Peak objects sorted by intensity (descending)
            
        Note:
            Only creates peaks above the threshold to reduce noise
        """
        peaks = []
        for mz, intensity in zip(mz_values, intensities):
            if intensity > threshold:
                peaks.append(Peak(mz=float(mz), intensity=float(intensity)))
        
        # Sort by intensity descending
        peaks.sort(key=lambda x: x.intensity, reverse=True)
        return peaks
    
    def _create_peaks_from_dataframe(self, df: pd.DataFrame) -> List[Peak]:
        """
        Create Peak objects from a pandas DataFrame.
        
        Args:
            df: DataFrame with 'mz' and 'intensity' columns
            
        Returns:
            List of Peak objects sorted by intensity (descending)
        """
        peaks = []
        for _, row in df.iterrows():
            peaks.append(Peak(mz=float(row['mz']), intensity=float(row['intensity'])))
        
        # Sort by intensity descending
        peaks.sort(key=lambda x: x.intensity, reverse=True)
        return peaks
    
    def _create_empty_peaks_result(self) -> Dict[str, Any]:
        """
        Create an empty peaks processing result structure.
        
        Used when no valid peak data is available, ensuring consistent
        output format for downstream processing.
        
        Returns:
            Dictionary with zero-filled arrays and metadata
        """
        max_peaks = self.config.max_peaks
        return {
            'peaks': np.zeros((max_peaks, 2), dtype=np.float32),
            'attention_mask': np.zeros(max_peaks, dtype=np.float32),
            'original_peak_count': 0,
            'processed_peak_count': 0,
            'max_mz': 0.0
        } 