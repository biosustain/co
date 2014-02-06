from __future__ import unicode_literals
from distutils.core import setup

setup(
    name='component-lib',
    version='0.1.0a1',
    packages=['colib', 'colib.utils'],
    url='',
    license='All Rights Reserved',
    author='Lars Schöning',
    author_email='lays@biosustain.dtu.dk',
    description='',
    test_suite='nose.collector',
    install_requires=[
        'biopython==1.62'
    ]
)
