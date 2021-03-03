
import glob
import os
import sys
import pytest
import subprocess

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
CIFS_DIR = os.path.join( os.path.join(TEST_DIR, os.path.pardir, 'cifs') )
CIF_FILES = glob.glob(os.path.join(CIFS_DIR, '*.cif'))

def run_cif2cell(args):
    return subprocess.check_output(['./binaries/cif2cell'] + args, stderr=subprocess.STDOUT).decode('utf8')

@pytest.mark.parametrize("cif_file", CIF_FILES)
def test_parse(cif_file):
    """Test running cif2cell on each CIF file in /cifs."""
    if sys.version_info < (3,0) and 'SiC.cif' in cif_file:
        pytest.skip(msg='skip test for files with unicode content under python 2.7.' +
           'see https://github.com/torbjornbjorkman/cif2cell/issues/7')

    result = run_cif2cell([cif_file])

    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result

def test_vasp():
    """Test VASP output."""
    cif_file = os.path.join(CIFS_DIR, "Si.cif")
    result = run_cif2cell(["-p", "vasp", "-f", cif_file])
    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result
    assert os.path.exists("POSCAR")

def test_castep():
    """Test CASTEP output."""
    cif_file = os.path.join(CIFS_DIR, "Si.cif")
    result = run_cif2cell(["-p", "castep", "-f", cif_file])

    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result
