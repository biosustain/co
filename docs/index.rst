.. co documentation master file, created by
   sphinx-quickstart on Mon Mar 17 16:34:51 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: co

Welcome to Co's documentation!
=================================

**Co** is a Python library for accessing and altering annotated DNA sequences. It keeps track of components and lifts
over feature annotations when a component is "mutated" by applying a series of mutations. With :mod:`co` you can
build new consensus sequences for cloned organisms and trace changes to features within a lineage.

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

- Future releases may include a *version tracking* system to track and propagate updated mutations due to e.g.
  better re-sequencing. Versioned components will maintain the relationship to the component's child component. Most
  likely, versions will be hashes of sequences and their mutations.
- An improved `non-strict` mode with better tolerance for overlapping mutations is planned.
- As this is a very early release of co, there is a long list of general improvements---they
  will be developed on demand.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

