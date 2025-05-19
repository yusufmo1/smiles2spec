"""
Feature extraction for molecules from SMILES strings.
"""

import numpy as np
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors, MACCSkeys, AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Avalon import pyAvalonTools

from .utils import handle_error, is_valid_smiles, logger

def extract_fingerprint(mol, fp_type, size=1024, radius=2, **kwargs):
    """Generic fingerprint extraction function."""
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
            return None
            
        fp = generator.GetFingerprint(mol)
        return np.array([int(b) for b in fp.ToBitString()])
    except Exception:
        return np.zeros(size, dtype=np.int32)

def extract_electronic_features(mol):
    """Extract electronic features from molecule."""
    features = {}
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
    
    features['electronic_features'] = np.array(list(electronic_features.values()))
    features['electronic_feature_names'] = list(electronic_features.keys())
    return features

def compute_molecular_features(smiles, config):
    """Compute all molecular features for a given SMILES string."""
    if not is_valid_smiles(smiles):
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        features = {}
        
        # Extract RDKit descriptors
        if config['extract_descriptors']:
            descriptor_names = [desc[0] for desc in Descriptors._descList]
            calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
            features['descriptors'] = np.array(calc.CalcDescriptors(mol))
            features['descriptor_names'] = descriptor_names
        
        # Extract fingerprints
        for fp_name, fp_config in config['fingerprints'].items():
            if fp_config.get('enabled', False):
                if fp_name == 'morgan':
                    for radius in fp_config['radii']:
                        fp_array = extract_fingerprint(
                            mol, 'morgan', size=fp_config['size'], radius=radius
                        )
                        features[f'morgan_fp_r{radius}'] = fp_array
                        if radius == fp_config['radii'][0]:
                            features['morgan_fingerprint'] = fp_array
                            # Add fingerprint bit names
                            features['morgan_fingerprint_names'] = [f'morgan_fingerprint_{i}' for i in range(len(fp_array))]
                elif fp_name == 'morgan_feature':
                    fp_array = extract_fingerprint(
                        mol, 'morgan_feature', size=fp_config['size'], radius=fp_config['radius']
                    )
                    features['morgan_feature_fp'] = fp_array
                    # Add fingerprint bit names
                    features['morgan_feature_fp_names'] = [f'morgan_feature_fp_{i}' for i in range(len(fp_array))]
                else:
                    fp_array = extract_fingerprint(
                        mol, fp_name, size=fp_config.get('size', 1024)
                    )
                    features[f'{fp_name}_fingerprint'] = fp_array
                    # Add fingerprint bit names
                    features[f'{fp_name}_fingerprint_names'] = [f'{fp_name}_fingerprint_{i}' for i in range(len(fp_array))]
        
        # Count features
        if config['count_features']['bond_counts']:
            bond_counts = defaultdict(int)
            for bond in mol.GetBonds():
                bond_type = str(bond.GetBondType())
                bond_counts[bond_type] += 1
            features['bond_counts'] = np.array(list(bond_counts.values()))
            features['bond_types'] = list(bond_counts.keys())
            
        if config['count_features']['atom_counts']:
            atom_counts = defaultdict(int)
            for atom in mol.GetAtoms():
                atom_counts[atom.GetSymbol()] += 1
            features['atom_counts'] = np.array(list(atom_counts.values()))
            features['atom_types'] = list(atom_counts.keys())
        
        # Extract electronic features
        if config['extract_electronic']:
            electronic_result = extract_electronic_features(mol)
            features.update(electronic_result)
            
        # Add SMILES to features
        features['smiles'] = smiles
        
        return features
    except Exception as e:
        error_result = handle_error(e, f"processing SMILES {smiles}", {"smiles": smiles})
        raise ValueError(f"Error processing SMILES: {error_result['error']}")