from __future__ import absolute_import
import glob
import os
import pytest
import subprocess

CIF_FILES = glob.glob('cifs/*.cif')

def run_cif2cell(args):
    return subprocess.check_output(['./binaries/cif2cell'] + args, stderr=subprocess.STDOUT).decode('utf8')

@pytest.mark.parametrize("cif_file", CIF_FILES)
def test_parse(cif_file):
    """Test running cif2cell on each CIF file in /cifs."""
    result = run_cif2cell([cif_file])

    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result

def test_vasp():
    """Test VASP output."""
    cif_file = "./cifs/Si.cif"
    result = run_cif2cell(["-p", "vasp", "-f", cif_file])
    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result
    assert os.path.exists("POSCAR")

def test_castep():
    """Test CASTEP output."""
    cif_file = "./cifs/Si.cif"
    result = run_cif2cell(["-p", "castep", "-f", cif_file])

    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result
