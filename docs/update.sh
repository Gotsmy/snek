#!/bin/bash
cd ..
pip uninstall snek -y
python setup.py install
cd docs
make html
# make latexpdf
