"""Feature extraction processor for molecules."""
import numpy as np
from collections import defaultdict
from typing import Dict, Any, List
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors, MACCSkeys, AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Avalon import pyAvalonTools

from ...config import settings
from ...core.models.molecule import MolecularFeatures, MoleculeInfo
from ...utils import logger, SMILESError

# Import utilities that will be migrated
from ...utils.chemistry import is_valid_smiles, smiles_to_png_base64, smiles_to_name

class FeatureProcessor:
    """
    Processor for extracting molecular features from SMILES strings.
    
    Extracts comprehensive molecular features including:
    - RDKit molecular descriptors (200+ properties)
    - Multiple types of molecular fingerprints
    - Electronic and charge distribution features
    - Atom and bond count statistics
    
    All features are computed using RDKit and optimized for
    mass spectrometry prediction tasks.
    
    Attributes:
        config: Feature extraction configuration from settings
        logger: Logger instance for error tracking
    """
    
    def __init__(self):
        self.config = settings.features
        self.logger = logger
    
    def extract_features(self, smiles: str) -> MolecularFeatures:
        """
        Extract comprehensive molecular features from a SMILES string.
        
        Computes all enabled feature types based on configuration settings.
        Features are used as input for mass spectrum prediction.
        
        Args:
            smiles: Valid SMILES string representation
            
        Returns:
            MolecularFeatures object containing:
                - descriptors: RDKit molecular descriptors array
                - fingerprints: Dictionary of various fingerprint types
                - electronic_features: Electronic property array
                - atom_counts: Dictionary of atom type counts
                - bond_counts: Dictionary of bond type counts
                
        Raises:
            SMILESError: If SMILES is invalid or feature extraction fails
        """
        if not is_valid_smiles(smiles):
            raise SMILESError(f"Invalid SMILES string: {smiles}")
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise SMILESError(f"Could not parse SMILES: {smiles}")
            
            # Extract descriptors
            descriptors = None
            descriptor_names = None
            if self.config.extract_descriptors:
                descriptors, descriptor_names = self._extract_descriptors(mol)
            
            # Extract fingerprints
            fingerprints = self._extract_fingerprints(mol)
            
            # Extract electronic features
            electronic_features = None
            if self.config.extract_electronic:
                electronic_features = self._extract_electronic_features(mol)
            
            # Extract counts
            atom_counts = {}
            bond_counts = {}
            if self.config.count_features.get('atom_counts', False):
                atom_counts = self._extract_atom_counts(mol)
            if self.config.count_features.get('bond_counts', False):
                bond_counts = self._extract_bond_counts(mol)
            
            return MolecularFeatures(
                smiles=smiles,
                descriptors=descriptors,
                descriptor_names=descriptor_names,
                fingerprints=fingerprints,
                electronic_features=electronic_features,
                atom_counts=atom_counts,
                bond_counts=bond_counts
            )
            
        except SMILESError:
            raise
        except Exception as e:
            self.logger.error(f"Error extracting features for SMILES {smiles}: {str(e)}")
            raise SMILESError(f"Feature extraction failed: {str(e)}", smiles)
    
    def extract_molecule_info(self, smiles: str) -> MoleculeInfo:
        """
        Extract basic molecular information for display purposes.
        
        Generates human-readable information about the molecule including
        chemical name, molecular weight, and structure visualization.
        
        Args:
            smiles: Valid SMILES string
            
        Returns:
            MoleculeInfo object containing:
                - chemical_name: IUPAC or common name (if available)
                - molecular_weight: Molecular weight in g/mol
                - exact_mass: Exact molecular mass
                - structure_png: Base64-encoded PNG of 2D structure
                
        Raises:
            SMILESError: If SMILES is invalid
        """
        if not is_valid_smiles(smiles):
            raise SMILESError(f"Invalid SMILES string: {smiles}")
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise SMILESError(f"Could not parse SMILES: {smiles}")
            
            # Get molecular weight and exact mass
            molecular_weight = Descriptors.MolWt(mol)
            exact_mass = Descriptors.ExactMolWt(mol)
            
            # Get chemical name
            chemical_name = smiles_to_name(smiles)
            
            # Get structure PNG
            structure_png = smiles_to_png_base64(smiles)
            
            return MoleculeInfo(
                smiles=smiles,
                chemical_name=chemical_name,
                molecular_weight=molecular_weight,
                exact_mass=exact_mass,
                structure_png=structure_png
            )
            
        except SMILESError:
            raise
        except Exception as e:
            self.logger.error(f"Error extracting molecule info for SMILES {smiles}: {str(e)}")
            raise SMILESError(f"Molecule info extraction failed: {str(e)}", smiles)
    
    def _extract_descriptors(self, mol) -> tuple[np.ndarray, List[str]]:
        """
        Extract all available RDKit molecular descriptors.
        
        Computes 200+ molecular properties including:
        - Physical properties (MW, LogP, TPSA, etc.)
        - Topological indices (Chi, Kappa indices)
        - Electronic properties (MR, atomic charges)
        - Fragment counts and functional groups
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Tuple of:
                - descriptors: NumPy array of descriptor values
                - descriptor_names: List of descriptor names in same order
                
        Note:
            Returns (None, None) if extraction fails
        """
        try:
            descriptor_names = [desc[0] for desc in Descriptors._descList]
            calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
            descriptors = np.array(calc.CalcDescriptors(mol))
            return descriptors, descriptor_names
        except Exception as e:
            self.logger.warning(f"Error extracting descriptors: {str(e)}")
            return None, None
    
    def _extract_fingerprints(self, mol) -> Dict[str, np.ndarray]:
        """
        Extract multiple types of molecular fingerprints.
        
        Generates various fingerprint representations based on configuration:
        - Morgan (circular) fingerprints with multiple radii
        - Morgan feature fingerprints (pharmacophore-based)
        - RDKit topological fingerprints
        - Avalon fingerprints
        - Pattern fingerprints
        - Layered fingerprints
        - MACCS keys
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Dictionary mapping fingerprint names to binary arrays
            
        Note:
            Only generates fingerprints enabled in configuration
        """
        fingerprints = {}
        
        for fp_name, fp_config in self.config.fingerprints.items():
            if not fp_config.get('enabled', False):
                continue
            
            try:
                if fp_name == 'morgan':
                    for radius in fp_config['radii']:
                        fp_array = self._extract_fingerprint(
                            mol, 'morgan', size=fp_config['size'], radius=radius
                        )
                        fingerprints[f'morgan_fp_r{radius}'] = fp_array
                        if radius == fp_config['radii'][0]:
                            fingerprints['morgan_fingerprint'] = fp_array
                elif fp_name == 'morgan_feature':
                    fp_array = self._extract_fingerprint(
                        mol, 'morgan_feature', 
                        size=fp_config['size'], 
                        radius=fp_config['radius']
                    )
                    fingerprints['morgan_feature_fp'] = fp_array
                else:
                    fp_array = self._extract_fingerprint(
                        mol, fp_name, size=fp_config.get('size', 1024)
                    )
                    fingerprints[f'{fp_name}_fingerprint'] = fp_array
            except Exception as e:
                self.logger.warning(f"Error extracting {fp_name} fingerprint: {str(e)}")
                continue
        
        return fingerprints
    
    def _extract_fingerprint(self, mol, fp_type: str, size: int = 1024, radius: int = 2) -> np.ndarray:
        """
        Extract a specific type of molecular fingerprint.
        
        Unified interface for generating different fingerprint types.
        Returns zero array on failure to ensure consistent dimensions.
        
        Args:
            mol: RDKit molecule object
            fp_type: Type of fingerprint to generate
            size: Bit vector size (not applicable to all types)
            radius: Radius for circular fingerprints (Morgan only)
            
        Returns:
            Binary fingerprint as NumPy array
            
        Supported Types:
            - morgan: ECFP-like circular fingerprints
            - morgan_feature: Feature-based circular fingerprints
            - topological: Topological torsion fingerprints
            - rdkit: RDKit path-based fingerprints
            - avalon: Avalon cheminformatics toolkit fingerprints
            - pattern: Substructure pattern fingerprints
            - layered: Layered fingerprints
            - maccs: MACCS structural keys (166 bits)
        """
        try:
            if fp_type == 'morgan':
                generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=size)
            elif fp_type == 'morgan_feature':
                generator = rdFingerprintGenerator.GetMorganGenerator(
                    radius=radius, 
                    fpSize=size,
                    atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
                )
            elif fp_type == 'topological':
                generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=size)
            elif fp_type == 'rdkit':
                return np.array([int(b) for b in Chem.RDKFingerprint(mol, fpSize=size).ToBitString()])
            elif fp_type == 'avalon':
                return np.array([int(b) for b in pyAvalonTools.GetAvalonFP(mol, nBits=size).ToBitString()])
            elif fp_type == 'pattern':
                return np.array([int(b) for b in Chem.PatternFingerprint(mol, fpSize=size).ToBitString()])
            elif fp_type == 'layered':
                return np.array([int(b) for b in Chem.LayeredFingerprint(mol, fpSize=size).ToBitString()])
            elif fp_type == 'maccs':
                return np.array(MACCSkeys.GenMACCSKeys(mol))
            else:
                return np.zeros(size, dtype=np.int32)
                
            fp = generator.GetFingerprint(mol)
            return np.array([int(b) for b in fp.ToBitString()])
        except Exception:
            return np.zeros(size, dtype=np.int32)
    
    def _extract_electronic_features(self, mol) -> np.ndarray:
        """
        Extract electronic and charge-related features.
        
        Computes features related to electronic structure and charge
        distribution, which are important for predicting fragmentation
        patterns in mass spectrometry.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            NumPy array of electronic features:
                - Partial charge statistics (min, max, mean, std, skew, kurtosis)
                - PEOE VSA descriptors (charge-weighted surface areas)
                - Crippen LogP and MR contributions
                
        Note:
            Uses Gasteiger charges as approximation of electron distribution.
            Returns zero values if calculation fails.
        """
        electronic_features = {}
        
        # Try to compute Gasteiger charges
        try:
            AllChem.ComputeGasteigerCharges(mol)
            charges = [atom.GetDoubleProp('_GasteigerCharge') if atom.HasProp('_GasteigerCharge') else 0.0 
                      for atom in mol.GetAtoms()]
            
            if charges:
                electronic_features['min_partial_charge'] = float(min(charges))
                electronic_features['max_partial_charge'] = float(max(charges))
                electronic_features['mean_partial_charge'] = float(np.mean(charges))
                electronic_features['charge_std'] = float(np.std(charges))
                
                if len(charges) > 1 and electronic_features['charge_std'] > 0:
                    electronic_features['charge_skew'] = float(np.mean(
                        [(c - electronic_features['mean_partial_charge'])**3 for c in charges]
                    ) / (electronic_features['charge_std']**3))
                    electronic_features['charge_kurtosis'] = float(np.mean(
                        [(c - electronic_features['mean_partial_charge'])**4 for c in charges]
                    ) / (electronic_features['charge_std']**4))
                else:
                    electronic_features['charge_skew'] = 0.0
                    electronic_features['charge_kurtosis'] = 0.0
            else:
                electronic_features.update({k: 0.0 for k in [
                    'min_partial_charge', 'max_partial_charge', 
                    'mean_partial_charge', 'charge_std', 
                    'charge_skew', 'charge_kurtosis'
                ]})
        except Exception:
            electronic_features.update({k: 0.0 for k in [
                'min_partial_charge', 'max_partial_charge', 
                'mean_partial_charge', 'charge_std', 
                'charge_skew', 'charge_kurtosis'
            ]})
        
        # Try to compute other electronic properties
        try:
            peoe_vsa = rdMolDescriptors.PEOE_VSA_(mol)
            electronic_features['PEOE_VSA1'] = peoe_vsa[1]
            electronic_features['PEOE_VSA2'] = peoe_vsa[2]
            electronic_features['PEOE_VSA3'] = peoe_vsa[3]
            
            crippen_contribs = rdMolDescriptors.GetCrippenContribs(mol)
            logp_values = [c[0] for c in crippen_contribs]
            mr_values = [c[1] for c in crippen_contribs]
            
            electronic_features['max_logp_contrib'] = float(max(logp_values)) if logp_values else 0.0
            electronic_features['min_logp_contrib'] = float(min(logp_values)) if logp_values else 0.0
            electronic_features['max_mr_contrib'] = float(max(mr_values)) if mr_values else 0.0
            electronic_features['min_mr_contrib'] = float(min(mr_values)) if mr_values else 0.0
        except Exception:
            electronic_features.update({k: 0.0 for k in [
                'PEOE_VSA1', 'PEOE_VSA2', 'PEOE_VSA3',
                'max_logp_contrib', 'min_logp_contrib',
                'max_mr_contrib', 'min_mr_contrib'
            ]})
        
        return np.array(list(electronic_features.values()))
    
    def _extract_atom_counts(self, mol) -> Dict[str, int]:
        """
        Count atoms by element type.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Dictionary mapping element symbols to counts
            Example: {'C': 6, 'H': 12, 'O': 1}
        """
        atom_counts = defaultdict(int)
        for atom in mol.GetAtoms():
            atom_counts[atom.GetSymbol()] += 1
        return dict(atom_counts)
    
    def _extract_bond_counts(self, mol) -> Dict[str, int]:
        """
        Count bonds by bond type.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Dictionary mapping bond types to counts
            Example: {'SINGLE': 5, 'DOUBLE': 2, 'AROMATIC': 6}
        """
        bond_counts = defaultdict(int)
        for bond in mol.GetBonds():
            bond_type = str(bond.GetBondType())
            bond_counts[bond_type] += 1
        return dict(bond_counts) 