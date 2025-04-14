
from rdkit import Chem
from rdkit.Chem import AllChem, MolStandardize
from rdkit.Chem.rdDistGeom import ETKDGv3
from rdkit.Chem.MolStandardize.rdMolStandardize import ChargeParent
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator


def parse_param_sets(param_sets_str):
    param_sets = []
    for param_set_str in param_sets_str.split(";"):
        params = {}
        for param in param_set_str.strip().split(","):
            key, value = param.strip().split("=")
            params[key.strip()] = value.strip().lower() == "true"
        # Generate a key for the parameter set
        key_parts = [f"+{k}" if v else f"-{k}" for k, v in params.items()]
        params["key"] = ",".join(key_parts)
        param_sets.append(params)
    return param_sets


def standardize_molecule(mol, *, stereo=True, isotope=True):
    # Clone the molecule to avoid modifying the original
    mol = Chem.Mol(mol)
    # Apply standardization procedures
    cleaner = MolStandardize.CleanupParameters()
    mol = MolStandardize.cleanup(mol, params=cleaner)
    # Optionally remove stereochemistry
    if not stereo:
        Chem.RemoveStereochemistry(mol)
    # Optionally remove isotope information
    if not isotope:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
    # Normalize the molecule
    normalizer = MolStandardize.normalize.Normalizer()
    mol = normalizer.normalize(mol)
    # Reionize the molecule
    reionizer = MolStandardize.reionize.Reionizer()
    mol = reionizer.reionize(mol)
    return mol


def generate_tautomers(mol, *, max_tautomers=1, min_score=0.0):
    enumerator = MolStandardize.tautomer.TautomerEnumerator()
    tautomers = enumerator.Enumerate(mol)
    # Score tautomers
    scored_tautomers = []
    for tautomer in tautomers:
        score = MolStandardize.tautomer.TautomerScoringFunctions.scoreTautomer(tautomer)
        if score >= min_score:
            scored_tautomers.append((tautomer, score))
    # Sort by score and limit the number
    scored_tautomers.sort(key=lambda x: x[1], reverse=True)
    return scored_tautomers[:max_tautomers]
