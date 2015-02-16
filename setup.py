# coding: utf-8
from __future__ import unicode_literals
import os
import io
from setuptools import setup, find_packages

def read(fname):
    return io.open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8').read()

setup(
    name='co',
    version='0.1.2',
    packages=find_packages(exclude=['*tests*']),
    url='http://co.readthedocs.org/en/latest/',
    license='MIT',
    author='Lars Schoening',
    author_email='lays@biosustain.dtu.dk',
    description='Python library for making and tracking mutated copies of DNA components',
    long_description=read('README.rst'),
    test_suite='nose.collector',
    install_requires=[
        'biopython>=1.63',
        'six>=1.5.2'
    ],
    setup_requires=[
        'nose>=1.1.2',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Other Environment',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    extras_require={
        'docs': ['Sphinx', 'sphinx-rtd-theme'],
    }
)
