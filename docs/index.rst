.. colib documentation master file, created by
   sphinx-quickstart on Mon Mar 17 16:34:51 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to colib's documentation!
=================================

``colib`` is a library for editing annotated DNA sequences.

The primary classes of importance in ``colib`` are :class:`Component`, :class:`Mutation` and :class:`HaploidOrganism`.
There is also a module for RESTful remote access to a *parts library* implemented in ``colib`` as well as sample
database models (these may eventually be moved into a separate project).

Contents:

.. toctree::
   :maxdepth: 2

   sequence


Road map
========

Current plans for future releases include a *version tracking* system of some sort that allows a mutated component to
be updated with a more accurate sequence while maintaining the relationship to the component's child component. Most
likely, versions will be hashes of the mutations that make up a component.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

