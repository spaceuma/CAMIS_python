# -*- coding: utf-8 -*-
from setuptools import setup

with open("README.md", 'r') as f:
    readme = f.read()
    
with open("LICENSE", 'r') as f:
    license = f.read()
    
setup(
      name='camis',
      version='0.1.0',
      description='CAMIS package for anisotropic cost modeling',
      long_description=readme,
      author='J. Ricardo Sanchez Ibanez',
      author_email='ricardosan@uma.es',
      url='https://github.com/spaceuma/CAMIS_python',
      license=license,
      packages=['camis']
      )
