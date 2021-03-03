# Copyright 2010 Torbjorn Bjorkman
# This file is part of cif2cell
#
# cif2cell is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cif2cell is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cif2cell.  If not, see <http://www.gnu.org/licenses/>.
#
from setuptools import setup, find_packages
from glob import glob

# Set up documentation
docfiles = ['docs/cif2cell.pdf']

# Get list of the cif example files
ciffiles = glob('cifs/*.cif')
periodiccifs = glob('cifs/periodic_table/*.cif')+['cifs/periodic_table/README']

setup(name='cif2cell',
      version='2.0.0a3',
      description='Construct a unit cell from CIF data',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      author='Torbjorn Bjorkman',
      author_email='torbjornb@gmail.com',
      url='http://cif2cell.sourceforge.net/',
      scripts=['binaries/cif2cell', 'binaries/vasp2cif'],
      python_requires=">=3.6",
      install_requires=[
          "PyCifRW==4.4",
      ],
      packages=find_packages(),
      data_files=[('lib/cif2cell', ['LICENSE']),
                  ('lib/cif2cell/sample_cifs', ciffiles),
                  ('lib/cif2cell/sample_cifs/periodic_table', periodiccifs),
                  ('lib/cif2cell/docs',docfiles)],
      license='GNU General Public License version 3',
      extras_require={
          'tests': [ 'pytest' ]
      },
      classifiers= [
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
    ],
)
