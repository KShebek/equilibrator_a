# Code to compute the Gibbs free energy of compounds given a SMILES
# Code is split into two areas and clearly marked
#   1. Code written by Kevin Shebek
#   2. Code taken from eQulibrator

from typing import Dict, List, Tuple, Optional

from pathlib import Path
import shutil
import os

import pandas as pd
import numpy as np
import quilt

from openbabel.pybel import readstring

from equilibrator_cache import Compound, CompoundMicrospecies, Q_
from equilibrator_cache.api import create_compound_cache_from_sqlite_file
from equilibrator_cache.compound_cache import PROTON_INCHI_KEY
from equilibrator_assets import chemaxon, thermodynamics

# The following modules were added from component_contribution/scripts/support
import group_decompose, molecule

# TODO: Use openbabel to reduce dependencies
from rdkit.Chem import AllChem

LOG10 = np.log(10.0)
CACHE_PATH = './cache'
SQLITE_PATH = Path('./cache/compounds.sqlite')
DEFAULT_QUILT_PKG = 'equilibrator/cache'

# Generate a compounds.sqlite file if it doesn't exist already
if not SQLITE_PATH.is_file():
    print('Local compounds.sqlite not found. Exporting from default equilibrator/cache quilt repo.')
    print('If another quilt repo is to be used, run generate_new_ccache with specified repo.')
    
    os.mkdir('./cache')
    quilt.install(DEFAULT_QUILT_PKG, force=True)
    quilt.export(DEFAULT_QUILT_PKG, './cache')

ccache = create_compound_cache_from_sqlite_file(SQLITE_PATH)
group_decomposer = group_decompose.GroupDecomposer()

###############################################################################
# KMS Code
def get_compound(mol_string: str, update_cache: bool = True, auto_commit: bool = False):
    """
    Gets a compound object from the compound cache, or generates one if not found.

    Parameters
    ----------
    mol_string : str
        A text description of the molecule(s) (SMILES or InChI). macOS is currently SMILES only.

    Returns
    -------   
    equilibrator_cache.Compound
        A Compound object that is used to calculate Gibbs free energy of formation and reactions.

    """
    # First check to see if compound is in ccache through partial InChI key match
    inchi_key = readstring('SMILES', mol_string).write("inchikey")
    cc_search = ccache.search_compound_by_inchi_key(inchi_key.split('-')[0])

    if cc_search:
        cpd = cc_search[0]
    else:
        cpd = _gen_compound(mol_string)
        # find next ID and add to sql database
        if update_cache == True:
            # Stage the compound for addition to the sqlite db
            ccache.session.add(cpd)
            # Flush executes the command queued up by session.add
            # temporary until .commit() is called
            ccache.session.flush()
            if auto_commit == True:
                ccache.session.commit()

    # assign an id if one wasn't already
    # need an id to calculate thermo properties
    if not cpd.id:
        cpd.id = -1

    return cpd

def _gen_compound(mol_string: str):
    """
    Generate an equilibrator_cache Compound object directly from a SMILES. 

    Parameters
    ----------
    mol_string : str
        A text description of the molecule(s) (SMILES).

    Returns
    -------   
    equilibrator_cache.Compound
        A Compound object that can be used to calculate Gibbs free energy of formation and reactions.

    """
    # TODO: This function currently only handles one mol_string input.

    mid_ph = 7
    # TODO: This simply gets around the inchi issue by passing whatever mol_string is to cxcalc. 
    # Should check OS and handle this accordingly.
    molecules = pd.DataFrame(
        data=[[-1, mol_string]],
        columns=["id", "inchi"]
    )

    # Calculate values to populate microspecies
    constants, pka_columns = chemaxon.get_dissociation_constants(
        molecules, '', num_acidic=20, num_basic=20, mid_ph=mid_ph
    )

    # Taken from equilibrator_assets/thermodynamics.py
    # Loops over constants dataframe and gets compound mappings
    min_ph = 0
    mid_ph = 7
    max_ph = 14
    compound_mappings = []
    for row in constants.itertuples(index=False):
        p_kas = [getattr(row, col) for col in pka_columns]
        p_kas = map(float, p_kas)
        p_kas = filter(lambda p_ka: min_ph < p_ka < max_ph, p_kas)
        dissociation_constants = sorted(p_kas, reverse=True)

        if pd.isnull(row.major_ms) or row.major_ms == "":
            compound_mappings.append(
                {
                    # "id": row.id,
                    "atom_bag": {},
                    "smiles": None,
                    "dissociation_constants": dissociation_constants,
                }
            )
        else:
            atom_bag = chemaxon.get_atom_bag("smi", row.major_ms)
            compound_mappings.append(
                {
                    # "id": row.id,
                    "atom_bag": atom_bag,
                    "smiles": row.major_ms,
                    "dissociation_constants": dissociation_constants,
                }
            ) 

    # Generate a compound with the compound mappings dictionary
    cpd = Compound(**compound_mappings[0])    

    # Specify cpd information not specified in compound_mappings
    major_ms = constants.iloc[0]['major_ms']
    mol = AllChem.MolFromSmiles(major_ms)
    cpd.inchi = AllChem.MolToInchi(mol)
    cpd.inchi_key = AllChem.MolToInchiKey(mol)

    # No magnesium data
    cpd.magnesium_dissociation_constants = []

    # Calculate microspecies and populate cpd with a list of CompoundMicrospecies generated with the microspecies dictionaries
    _, microspecies = _get_microspecies_data(cpd.id, major_ms, cpd.dissociation_constants, cpd.atom_bag)
    cpd.microspecies = [CompoundMicrospecies(**ind_ms) for ind_ms in microspecies]    

    # Decompose the compounds into the group vectors
    mol = molecule.Molecule.FromSmiles(major_ms)
    decomposition = group_decomposer.Decompose(mol, ignore_protonations=False, raise_exception=True)
    cpd.group_vector = list(decomposition.AsVector())
      
    return cpd

def generate_new_ccache(
    package: str = DEFAULT_QUILT_PKG,
    hash: Optional[str] = None,
    tag: Optional[str] = None,
    version: Optional[str] = None,
    force: bool = True,
    ):
    """
    Generate the new compound.sqlite from a specified quilt package

    Parameters
    ----------
    package : str, optional
        The quilt package used to initialize the compound cache
        (Default value = `equilibrator/cache`)
    hash : str, optional
        quilt hash (Default value = None)
    version : str, optional
        quilt version (Default value = None)
    tag : str, optional
        quilt tag (Default value = None)
    force : bool, optional
        Re-download the quilt data if a newer version exists
        (Default value = `True`).

    """

    print(f"Attempting to install {package} from quilt.")
    quilt.install(
        package=package, hash=hash, tag=tag, version=version, force=force
        )
    print('Removing old cache and generate new one')
    shutil.rmtree(CACHE_PATH, ignore_errors=True)
    os.mkdir('./cache')
    quilt.export(package, './cache')
    ccache = create_compound_cache_from_sqlite_file(SQLITE_PATH)   

def _get_ccache():
    """
    Returns the ccache object for direct use. Use with caution if you aren't very familiar with how the
    eQuilibrator sqlite database works!

    Returns
    ----------
    ccache : equilibrator_cache.compound_cache.CompoundCache
        The compound cache generated from compounds.sqlite
    """
    return ccache

# End KMS Code
###############################################################################

###############################################################################
# Noor Code
def _get_microspecies_data(
    cpd_id: int,
    major_ms: int,
    dissociation_constants: List[float],
    atom_bag: Dict[str, int],
    mid_ph: float = 7.0,
) -> Tuple[dict, List[dict]]:
    """
    Calculate the microspecies information for a compound (if possible).

    Returns
    -------
    tuple
        dict
            A mapping for updating a compound with atom bag and dissociation
            constants.
        list
            A list of microspecies mappings for that compound.

    """
    # Compounds for which the major microspecies calculation failed are skipped.
    if pd.isnull(major_ms) or major_ms == "":
        return (
            {
                "id": id,
                "atom_bag": {},
                "smiles": None,
                "dissociation_constants": dissociation_constants,
            },
            [],
        )

    major_ms_num_protons = atom_bag.get("H", 0)
    num_protons = sum(count * chemaxon.SYMBOL_TO_ATOMIC_NUMBER[elem]
                      for elem, count in atom_bag.items() if elem != "e-")
    major_ms_charge = num_protons - atom_bag.get("e-", 0)

    num_species = len(dissociation_constants) + 1
    # Find the index of the major microspecies, by counting how many pKas there
    # are in the range between the given pH and the maximum (typically, 7 - 14).

    # KMS .any()
    if not dissociation_constants:
        major_ms_index = 0
        num_protons = [major_ms_num_protons]
        charges = [major_ms_charge]
    else:
        major_ms_index = sum(
            (1 for p_ka in dissociation_constants if p_ka > mid_ph)
        )
        num_protons = [
            i - major_ms_index + major_ms_num_protons
            for i in range(num_species)
        ]
        charges = [
            i - major_ms_index + major_ms_charge for i in range(num_species)
        ]

    microspecies = []
    for i, (z, nH) in enumerate(zip(charges, num_protons)):
        is_major = False

        if i == major_ms_index:
            ddg_over_rt = 0.0
            is_major = True
        elif i < major_ms_index:
            ddg_over_rt = sum(dissociation_constants[i:major_ms_index]) * LOG10
        elif i > major_ms_index:
            ddg_over_rt = -sum(dissociation_constants[major_ms_index:i]) * LOG10
        else:
            raise IndexError("Major microspecies index mismatch.")
        microspecies.append(
            {
                "compound_id": cpd_id,
                "charge": z,
                "number_protons": nH,
                "ddg_over_rt": ddg_over_rt,
                "is_major": is_major,
                "number_magnesiums": 0
            }
        )
    return (
        {
            "id": cpd_id,
            "atom_bag": atom_bag,
            "smiles": major_ms,
            "dissociation_constants": dissociation_constants,
        },
        microspecies,
    )

def _populate_microspecies(session, compound: Compound, mid_ph: float = 7.0) -> None:
    """
    Calculate dissociation constants and create microspecies for a single Compound.

    Parameters
    ----------
    session : sqlalchemy.orm.session.Session
        An active session in order to communicate with a SQL database.
    compound: Compound
        The compound object to determine the microspecies.
    mid_ph : float
        The pH for which the major microspecies is calculated
        (Default value = 7.0).

    """

    # We only create microspecies for compounds that have dissociation_constants
    # (although it could also be an empty list) and an atom_bag (which we need
    # in order to determine the nH and z of the major microspecies).

    microspecies_mappings = _create_microspecies_mappings(compound, mid_ph)
    session.bulk_insert_mappings(CompoundMicrospecies, microspecies_mappings)


def _create_microspecies_mappings(
    compound: Compound, mid_ph: float = 7.0
) -> List[dict]:
    """Create the mappings for the microspecies of a Compound.

    Parameters
    ----------
    compound : Compound
        A Compound object, where the atom_bag and dissociation_constants must
        not be None.
    mid_ph : float
        The pH for which the major microspecies is calculated
        (Default value = 7.0).

    Returns
    -------
    list
        A list of mappings for creating the entries in the compound_microspecies
        table.

    """
    # We add an exception for H+ (and put z = nH = 0) in order to
    # eliminate its effect of the Legendre transform.
    if compound.inchi_key == PROTON_INCHI_KEY:
        return [
            {
                "compound_id": compound.id,
                "charge": 0,
                "number_protons": 0,
                "number_magnesiums": 0,
                "ddg_over_rt": 0.0,
                "is_major": True,
            }
        ]

    # Find the index of the major microspecies, by counting how many pKas there
    # are in the range between the given pH and the maximum (typically, 7 - 14).
    # Then make a list of the nH and charge values for all the microspecies
    if not compound.dissociation_constants:
        num_species = 1
        major_ms_index = 0
    else:
        num_species = len(compound.dissociation_constants) + 1
        major_ms_index = sum(
            (1 for p_ka in compound.dissociation_constants if p_ka > mid_ph)
        )

    major_ms_num_protons = compound.atom_bag.get("H", 0)
    major_ms_charge = compound.net_charge

    microspecies_mappings = dict()
    for i in range(num_species):
        charge = i - major_ms_index + major_ms_charge
        num_protons = i - major_ms_index + major_ms_num_protons

        if i == major_ms_index:
            ddg_over_rt = 0.0
        elif i < major_ms_index:
            ddg_over_rt = (
                sum(compound.dissociation_constants[i:major_ms_index]) * LOG10
            )
        elif i > major_ms_index:
            ddg_over_rt = (
                -sum(compound.dissociation_constants[major_ms_index:i]) * LOG10
            )
        else:
            raise IndexError("Major microspecies index mismatch.")

        microspecies_mappings[(num_protons, 0)] = {
            "compound_id": compound.id,
            "charge": charge,
            "number_protons": num_protons,
            "number_magnesiums": 0,
            "ddg_over_rt": ddg_over_rt,
            "is_major": i == major_ms_index,
        }

    return sorted(
        microspecies_mappings.values(),
        key=lambda x: (x["number_magnesiums"], x["number_protons"]),
    )

# End Noor Code
###############################################################################