# coding: utf-8
from __future__ import unicode_literals
from setuptools import setup, find_packages

setup(
    name='colib',
    version='0.1.0',
    packages=find_packages(exclude=['*tests*']),
    url='',
    license='MIT',
    author='Lars Schöning',
    author_email='lays@biosustain.dtu.dk',
    description='',
    test_suite='nose.collector',
    install_requires=[
        'Flask-Presst>=0.2.2',
        'biopython>=1.63'
    ],
    extras_require={
        'docs': ['Sphinx', 'sphinx-rtd-theme'],
    }
)
