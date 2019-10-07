import unittest
import subprocess
import glob
import os


class TestCif2Cell(unittest.TestCase):
    def test_cif(self):
        for path in glob.glob('cifs/*.cif'):
            p = subprocess.Popen(['./binaries/cif2cell', path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout_data, stderr_data = p.communicate()
            print("OUT")
            print(stdout_data.decode())

            assert not "***Warning: Space group operation check failed" in stdout_data.decode()
            assert not "Error" in stdout_data.decode()

    def test_vasp(self):
        # for path in glob.glob('cifs/*.cif'):
        p = subprocess.Popen(['./binaries/cif2cell', "-p", "vasp", "-f", 'cifs/Si.cif'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_data, stderr_data = p.communicate()
        print("OUT")
        print(stdout_data.decode())

        assert not "***Warning: Space group operation check failed" in stdout_data.decode()
        assert not "Error" in stdout_data.decode()
        assert os.path.exists("POSCAR")

    def test_castep(self):
        # for path in glob.glob('cifs/*.cif'):
        p = subprocess.Popen(['./binaries/cif2cell', "-p", "castep", "-f", 'cifs/Si.cif'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_data, stderr_data = p.communicate()
        print("OUT")
        print(stdout_data.decode())

        assert not "***Warning: Space group operation check failed" in stdout_data.decode()
        assert not "Error" in stdout_data.decode()