# coding: utf-8
from __future__ import unicode_literals
from setuptools import setup

setup(
    name='component-lib',
    version='0.1.0a1',
    packages=['colib', 'colib.utils'],
    url='',
    license='MIT',
    author='Lars Schöning',
    author_email='lays@biosustain.dtu.dk',
    description='',
    test_suite='nose.collector',
    install_requires=[
        'Flask-Presst>=0.2.2',
        'biopython>=1.62',
        'blinker>=1.3',
    ]
)
