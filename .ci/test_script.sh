#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# try cifs
cd cifs
for cif in *.cif; do
    cif2cell $cif
done

# try periodic table
cd periodic_table
for cif in *.cif; do
    cif2cell $cif
done
