from typing import Dict, List, Tuple, Optional

from pathlib import Path
import shutil
import os

import pandas as pd
import numpy as np
import quilt

from openbabel.pybel import readstring

from equilibrator_cache.api import create_compound_cache_from_sqlite_file
from equilibrator_assets.generate_compound import create_compound

LOG10 = np.log(10.0)
CACHE_PATH = './cache'
SQLITE_PATH = Path('./cache/compounds.sqlite')
DEFAULT_QUILT_PKG = 'equilibrator/cache'

# Generate a compounds.sqlite file if it doesn't exist already
if not SQLITE_PATH.is_file():
    print('Local compounds.sqlite not found. Exporting from default equilibrator/cache quilt repo.')
    print('If another quilt repo is to be used, run generate_new_ccache with specified quilt location.')
    
    os.mkdir('./cache')
    quilt.install(DEFAULT_QUILT_PKG, force=True)
    quilt.export(DEFAULT_QUILT_PKG, './cache')

ccache = create_compound_cache_from_sqlite_file(SQLITE_PATH)

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
    cc_search = ccache.search_compound_by_inchi_key(inchi_key.rsplit('-', 1)[0])

    if cc_search:
        cpd = cc_search[0]
    else:
        cpd = create_compound(mol_string)
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
    Returns the ccache object for direct use. 

    Returns
    ----------
    ccache : equilibrator_cache.compound_cache.CompoundCache
        The compound cache generated from compounds.sqlite
    """
    return ccache