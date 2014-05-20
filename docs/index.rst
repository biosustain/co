.. colib documentation master file, created by
   sphinx-quickstart on Mon Mar 17 16:34:51 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: colib

Welcome to colib's documentation!
=================================

``colib`` is a library for altering DNA sequence records in BioPython. Using mutations, the library lets you create ne
 sequence records and lift over features from one sequence to another.

The primary classes of importance in ``colib`` are :class:`Component`, :class:`mutations.Mutation`
and :class:`HaploidOrganism`. These can encode lineages of annotated DNA components and strains.

Contents:
---------

.. toctree::
   :maxdepth: 2

   quickstart
   translation
   mutations
   component


Road map
========

Current plans for future releases include a *version tracking* system of some sort that allows a mutated component to
be updated with a more accurate sequence while maintaining the relationship to the component's child component. Most
likely, versions will be hashes of the mutations that make up a component.

Also planned is a *context* interface that can be implemented to automatically save, load and look up parts in the
context of a storage backend.

A improved `non-strict` mode with better tolerance for overlapping mutations is also being considered.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

