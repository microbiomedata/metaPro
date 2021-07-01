'''
Helps to resolve the dependency issues among packages.! In order to import from sibling packages.
You don't have to add every sibling package to "sys.path", instead install your project itself in editable state!
      - Create a virtual environment if needed.
      - add a setup file as below
      - Install your project in editable state!
        $ pip install -e .

Note: Run this script whenever you add to the directory structure.
'''
from setuptools import setup, find_packages
setup(name = 'proteomics_workflow',
      version= '0.1',
      packages = find_packages())