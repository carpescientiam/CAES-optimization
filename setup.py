#! /usr/bin/env python
from setuptools import setup

setup(
    name='caes',
    version='0.0.2',
    author='Fahim Sadat',
    author_email='fahim.sadat@hs-flensburg.de',
    description='Models for compressed air energy storage (CAES)',
    package_dir={'caes': 'caes'},
    install_requires=['scipy == 1.10.0', 'numpy >= 1.14.0']
)
