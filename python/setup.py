#!/usr/bin/env python
import sys

if 'develop' in sys.argv:
    # use setuptools for develop, but nothing else
    from setuptools import setup
else:
    from distutils.core import setup

with open('README.rst') as file:
    long_description = file.read()

setup(name='spectrum_extraction',
      version='0.1',
      description='Tools to extract spectra from data cubes given a mask',
      long_description=long_description,
      author='Adam Ginsburg',
      author_email='adam.g.ginsburg@gmail.com',
      url='https://github.com/BGPS/distance_omnibus2',
      packages=['spectrum_extraction','pca_distance','nirex'],
      scripts=['scripts/*py'],
     )
