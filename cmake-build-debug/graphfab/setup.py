#! /usr/bin/env python

from setuptools import setup

setup(name='sbnw',
      version='1.3.27',
      description='A package for automatically laying out SBML models',
      author='J Kyle Medley',
      url='https://github.com/0u812/sbnw',
      packages=['sbnw'],
      package_data={'sbnw': ['*.so*','*.dylib*','lib*']})
