# coding: utf-8
from __future__ import unicode_literals
from setuptools import setup, find_packages

setup(
    name='co',
    version='0.1.0',
    packages=find_packages(exclude=['*tests*']),
    url='',
    license='MIT',
    author='Lars Schoening',
    author_email='lays@biosustain.dtu.dk',
    description='',
    test_suite='nose.collector',
    install_requires=[
        'biopython>=1.63'
    ],
    extras_require={
        'docs': ['Sphinx', 'sphinx-rtd-theme'],
    }
)
