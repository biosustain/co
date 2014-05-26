.. colib documentation master file, created by
   sphinx-quickstart on Mon Mar 17 16:34:51 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: colib

Welcome to colib's documentation!
=================================

COLIB is a Python library for accessing and altering annotated DNA sequences. It keeps track of components and lifts
over feature annotations when a component is "mutated" by applying a series of mutations. With :mod:`colib` you can
build new consensus sequences for cloned organisms and follow changes to feature annotations within a lineage.

Contents:
---------

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   mutations
   component
   organisms
   translation

Road map
========

Current plans for future releases include a *version tracking* system of some sort that allows a mutated component to
be updated with a more accurate sequence while maintaining the relationship to the component's child component. Most
likely, versions will be hashes of the mutations that make up a component.

An improved `non-strict` mode with better tolerance for overlapping mutations is also being considered.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

