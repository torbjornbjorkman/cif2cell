"""Tests for cif2cell."""
from pathlib import Path
import subprocess
import platform
import os
import pytest

TEST_DIR = Path(__file__).resolve().parent
CIFS_DIR = TEST_DIR.parent / 'cifs'
CIF_FILES = CIFS_DIR.glob('*.cif')
CIF2CELL_SCRIPT = TEST_DIR.parent / 'binaries' / 'cif2cell'

def run_cif2cell(args):
    if platform.system() == 'Windows':
        # windows does not understand the shebang
        cmd = ["python.exe", CIF2CELL_SCRIPT] + args
    else:
        cmd = [ CIF2CELL_SCRIPT ] + args

    try:
        # Setting the PYTHONIOENCODING was necessary for windows runners on GitHub action
        result = subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT,
            env=dict(os.environ, PYTHONIOENCODING='utf8'),
            encoding='utf8')
        return result
    except subprocess.CalledProcessError as exc:
        raise EnvironmentError(result) from exc


@pytest.mark.parametrize("cif_file", CIF_FILES)
def test_parse(cif_file):
    """Test running cif2cell on each CIF file in /cifs."""
    result = run_cif2cell([cif_file])

    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result

def test_vasp():
    """Test VASP output."""
    cif_file = CIFS_DIR / "Si.cif"
    result = run_cif2cell(["-p", "vasp", "-f", cif_file])
    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result
    assert Path("POSCAR").exists()

def test_castep():
    """Test CASTEP output."""
    cif_file = CIFS_DIR / "Si.cif"
    result = run_cif2cell(["-p", "castep", "-f", cif_file])

    assert not "***Warning: Space group operation check failed" in result
    assert not "Error" in result
