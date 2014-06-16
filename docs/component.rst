
=====================
Components & Features
=====================

.. module:: co

The :class:`Component` class is very similar to :class:`Bio.SeqRecord.SeqRecord`, but does not currently
sub-class it---mainly because the ``features`` property is implemented differently.


.. automodule:: co
   :members: Component
   :undoc-members:


The :class:`Feature` class inherits from :class:`Bio.SeqFeature.SeqFeature` but stores some additional information.
Proceed with caution when using the two types interchangeably.



.. automodule:: co
   :members: Feature, FeatureSet, ComponentFeatureSet

