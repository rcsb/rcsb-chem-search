# Correct RDKit imports
import base64
import operator
import time
from dataclasses import dataclass, field
from collections.abc import Callable, Iterator, Generator
from typing import Self

import numpy as np
from loguru import logger
from rdkit import Chem
from rdkit.Chem import Mol, AssignCIPLabels
from rdkit.Chem.rdDistGeom import ETKDGv3, EmbedParameters,EmbedMultipleConfs
from rdkit.Chem.rdMolTransforms import CanonicalizeMol # canonicalizes conformers
import rdkit.rdBase as rdBase
from rdkit.Chem.MolStandardize.rdMolStandardize import (
    TautomerEnumerator,
    CleanupParameters,
    Normalizer,
    Reionizer,
    TautomerEnumeratorResult,
    TautomerEnumeratorStatus,
)
from rdkit.Chem import rdmolops

from .core.errors import AnyChemError, MoleculeBuildError, MoleculeOperationError
from .models import Structure, Tautomer, OpStatus, MoleculeFactory
from .core.rdkit_utils import RdkitLogging

_state_map: dict[int, OpStatus] = {
    TautomerEnumeratorStatus.Completed: OpStatus.COMPLETED,
    TautomerEnumeratorStatus.MaxTransformsReached: OpStatus.HALTED,
    TautomerEnumeratorStatus.Failed: OpStatus.FAILED,
}
_max_tautomer_score = 1_000_000

rdmolops.SetUseLegacyStereoPerception(False)
rdmolops.SetUseLegacyStereoPerception(False)
# See: DativeBondsToHaptic and HapticBondsToDative
# See: FindPotentialStereoBonds
# See: GetFormalCharge
# See: Kekulize and KekulizeIfPossible
# See: RemoveStereochemistry
# See: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.SanitizeFlags


#@Clz.describe("{<smiles>} score={score}")
@dataclass(slots=True, frozen=True)
class EnumeratedTautomer:
    smiles: str
    mol: Mol
    score: int


@dataclass(slots=True, frozen=True)
class EnumeratedTautomerResults:
    state: OpStatus
    tautomers: list[EnumeratedTautomer]
    n_enumerated: int
    n_passed: int
    n_limit: int
    n_modified_atoms: int
    n_modified_bonds: int
    seconds_taken: float
    errors: list[Exception]

    @property
    def n_contained(self) -> int:
        return len(self.tautomers)

    def __getitem__(self, index: int) -> EnumeratedTautomer:
        return self.tautomers[index]


@dataclass(frozen=True)
class TautomerEnumeratorFactory:

    def create(
        self,
        reassign_stereo: bool = False,
        max_initial: int | None = None,
        max_transforms: int = 1000
    ) -> TautomerEnumerator:
        # https://github.com/rdkit/rdkit/blob/2e5f7ce80c4c9ae11c8573593f4a08c34c3f5166/Code/GraphMol/MolStandardize/Tautomer.h#L196
        enumerator = TautomerEnumerator()
        enumerator.SetMaxTransforms(max_transforms)
        if max_initial is not None:
            enumerator.setMaxTautomers(max_initial)
        # TautomerEnumerator tries to handle stereochemistry correctly.
        # The code loops over all atoms and removes iff the atom is SP2 OR removeSp3Stereo is set.
        # https://github.com/rdkit/rdkit/blob/2e5f7ce80c4c9ae11c8573593f4a08c34c3f5166/Code/GraphMol/MolStandardize/Tautomer.cpp#L169
        # It then does something similar for all bonds:
        # https://github.com/rdkit/rdkit/blob/2e5f7ce80c4c9ae11c8573593f4a08c34c3f5166/Code/GraphMol/MolStandardize/Tautomer.cpp#L198
        # So, these options (removeSp3Stereo and removeBondStereo) really do remove ALL stereochemistry.
        enumerator.SetRemoveBondStereo(reassign_stereo)
        enumerator.SetRemoveSp3Stereo(reassign_stereo)
        # The reassignStereo option uses AssignStereochemistry, which is inaccurate but fast.
        # CIPLabeler is more accurate but slower.
        # https://www.rdkit.org/docs/source/rdkit.Chem.rdCIPLabeler.html
        enumerator.SetReassignStereo(False)
        #enumerator.setRemoveIsotopicHs(False)
        return enumerator


@dataclass(slots=True, frozen=True, eq=False)
class TautomerGenerator:
    molecule_factory: MoleculeFactory
    enumerator: TautomerEnumerator = field(default_factory=lambda: TautomerEnumerator(), init=False)
    min_score: int = 0
    max_tautomers: int = 100
    max_cip_iterations: int = 2E6  # 0.1 ms; only likely for peptides.

    def __post_init__(self) -> None:
        # https://www.rdkit.org/docs/cppapi/MolStandardize_2Tautomer_8h_source.html
        # use maxTransforms rather than maxTautomers
        # the former scales more linearly w.r.t. the number of tautomeric centers

    def apply(self, source_mol: Mol, source_smiles: str) -> EnumeratedTautomerResults:
        t0 = time.monotonic()
        with RdkitErrors.wrap(op_name="enumerate tautomers"):
            result: TautomerEnumeratorResult = self.enumerator.Enumerate(source_mol)
        state: OpStatus = _state_map[result.status]
        passed: list[EnumeratedTautomer] = [EnumeratedTautomer(source_smiles, source_mol, _max_tautomer_score)]
        errors: list[MoleculeOperationError] = []
        if state.can_contain_results:
            for tautomer_smiles, tautomer_mol in result.smilesTautomerMap.items():
                match (result := self._add(tautomer_smiles, tautomer_mol)):
                    case EnumeratedTautomer():
                        passed.append(result)
                    case MoleculeOperationError():
                        if result.score >= self.min_score:
                            errors.append(result)
            passed.sort(key=lambda t: t.score, reverse=True)  # sort highest score first
        final_hits: list[EnumeratedTautomer] = passed[:self.max_tautomers]
        return EnumeratedTautomerResults(
            state=state,
            tautomers=final_hits,
            n_enumerated=len(result.smilesTautomerMap),
            n_passed=len(passed),
            n_limit=self.max_tautomers,
            n_modified_atoms=len(result.modifiedAtoms),
            n_modified_bonds=len(result.modifiedBonds),
            seconds_taken=time.monotonic() - t0,
            errors=errors,
        )

    def _add(self, tautomer_mol, tautomer_smiles: str) -> EnumeratedTautomer | MoleculeOpError:
        # self.enumerator.tautomerScoreVersion
        # "inspired by" https://doi.org/10.1007/s10822-010-9346-4
        # https://github.com/rdkit/rdkit/blob/9e2b3f233eafdb30fe51a57b283438a57a537719/Code/GraphMol/MolStandardize/Tautomer.cpp#L35
        try:
            with RdkitErrors.wrap(op_name="score tautomer"):
                score = self.enumerator.ScoreTautomer(tautomer_mol)
        except MoleculeOperationError as e:
            return e
        return EnumeratedTautomer(tautomer_smiles, tautomer_mol, score)


@dataclass(slots=True, frozen=True, eq=False)
class Standardizer:
    normalizer: Normalizer = field(default_factory=lambda: Normalizer(), init=False)
    reionizer: Reionizer = field(default_factory=lambda: Reionizer(), init=False)
    cleanup: CleanupParameters = field(default_factory=lambda: CleanupParameters(), init=False)

    def standardize(self, mol: Mol) -> Mol:
        with RdkitErrors.wrap(op_name="normalize"):
            mol = self.normalizer.normalize(mol)
        with RdkitErrors.wrap(op_name="re-ionize"):
            mol = self.reionizer.reionize(mol)
        with RdkitErrors.wrap(op_name="re-ionize"):
            pass  # TODO mol = self.**(mol)
        return mol


@dataclass(slots=True, frozen=True)
class LigandResults:
    ligand: Structure
    tautomers: list[Tautomer]


@dataclass(slots=True, frozen=True, eq=False)
class Pipeline:
    molecule_factory: MoleculeFactory
    tautomer_generator: TautomerGenerator
    standardizer: Standardizer
    fingerprints: set[str]
    fingerprinter: Callable[[str, Mol], np.array]

    def create_dataset(self, source: Iterator[tuple[str, str]]) -> None:
        for rcsb_id, smiles in source:
            ligand = self._get_std(rcsb_id, smiles)
            tautomers = [
                self.molecule_factory.of(t.mol,key=f"{rcsb_id}.{i}")
                for i, t in enumerate(self._tautomers(ligand))
            ]
            ligand = self.molecule_factory.of(smiles, key=rcsb_id)
            results = self.process(rcsb_id, smiles)

    def _get_std(self,rcsb_id,smiles: str)->Structure:
        ligand = self.molecule_factory.of(smiles)
        mol = self.standardizer.standardize(ligand.mol)
        return self.molecule_factory.of(mol, key=rcsb_id)

    def encode_fingerprint(self, fp_array: np.array) -> str:
        # Convert to bytes and encode
        fp_bytes = np.packbits(fp_array.astype(np.int16))
        encoded_fp = base64.b64encode(fp_bytes).decode("utf-8")
        return encoded_fp

    def process(self, ligand:Structure) -> list[Tautomer]:
        pass

    def _tautomers(self, ligand: Structure) -> Generator[Tautomer]:

        def log(level: str, msg: str, *args, **kwargs) -> None:
            logger.log(level, "[{rcsb_id}/{inchikey}]>> " + msg, ligand.key, ligand.inchikey, *args, **kwargs)

        results: EnumeratedTautomerResults = self.tautomer_generator.apply(ligand.mol)
        log(
            "DEBUG",
            "{u} used of {e} enumerated ({p} passed)",
            u=results.n_contained,
            e=results.n_enumerated,
            p=results.n_passed,
        )
        n_errors: int = len(results.errors)
        log("DEBUG", "took {s:.2f} s", s=results.seconds_taken)
        if Log.is_debug_enabled() and n_errors > 0:
            log("DEBUG", "encountered {n} errors", n=n_errors)
            for i, e in enumerate(results.errors):
                log("DEBUG", "error {i}/{n}: {e}", i, n=n_errors, e=e, exc_info=True)
        if Log.is_trace_enabled():
            for i, tautomer in results.tautomers:
                logger.trace("tautomer {i}/{n}: {t}", i=i, n=results.n_contained, t=repr(tautomer))
        # convert EnumeratedTautomer to Tautomer
        for i, tautomer in enumerate(results.tautomers):
            yield Tautomer(
                **self.molecule_factory.of(tautomer.mol).as_dict,
                score=tautomer.score,
            )
