
==========
Quickstart
==========

Simple mutation
---------------

.. module:: colib

.. without features

To illustrate what `colib` is designed for, let's begin with a hello world example:

.. code-block:: python

    >>> from colib import Component
    >>> from colib.mutations import *
    >>> hi_x = Component('Hello X!')
    >>> hi_x.seq
    Seq('Hello X!', Alphabet())
    >>> hi_world = hi_x.mutate([Mutation(6, 1, 'world')])
    >>> hi_world.seq
    Seq('Hello world!', Alphabet())



- components
- mutations
- features


.. module:: mutation

Mutation types that inherit from :class:`Mutation` are :class:`INS`,  :class:`DEL`, :class:`SUB` for substitutions, and
:class:`SNP` for SNPs.

.. module:: colib
.. features

Features & feature inheritance
------------------------------

.. code-block:: python

    >>> from Bio.SeqFeature import *
    >>>
    >>> slogan = Component('Colib is for DNA components')
    >>> slogan.features.add(FeatureLocation(0, 5), type='library')
    Feature(FeatureLocation(ExactPosition(0), ExactPosition(5)), type='library')
    >>> slogan.features.add(FeatureLocation(13, 16), id='DNA')
    Feature(FeatureLocation(ExactPosition(13), ExactPosition(16)), id='DNA')
    >>>
    >>> slogan.features.add(FeatureLocation(17, 28)).seq
    Seq('components', Alphabet())
    >>>
    >>> [f.seq for f in slogan.features]
    [Seq('Colib', Alphabet()), Seq('DNA', Alphabet()), Seq('components', Alphabet())]
    >>>


FIXME Colib should not become Colibary, but Co[lib] could become Co[library]?


`colib` also translates feature annotations when mutating:


.. code-block:: python

    >>> new_slogan = slogan.mutate([DEL(13, 3), INS(5, 'rary')])
    >>> new_slogan.seq
    Seq('Colibrary is for  components', Alphabet())
    >>> new_slogan.features
    ComponentFeatureSet([Feature(FeatureLocation(ExactPosition(0), ExactPosition(9)), type='library'),
                         Feature(FeatureLocation(ExactPosition(19), ExactPosition(29)))])
    >>>
    >>>
    >>> [f.seq for f in new_slogan.features]
    [Seq('Colibrary', Alphabet()), Seq('components', Alphabet())]



Optimization behind the scenes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Feature annotations that are inherited from another component are not copied over
in memory --- instead they are looked up on the fly. Only added and removed features are stored. A feature is
considered changed when its sequence is affected in any way. When a feature is changed, the old feature is removed and
the new feature is added.

- On-the-fly coordinate translation of unchanged features is done using :class:`translation.TranslationTable`, inspired
  by the UCSC's liftOver format.
- Feature locations are indexed using :class:`interval.IntervalTree`, currently implemented as a BST.


Feature diffs
^^^^^^^^^^^^^



.. code-block:: python

    >>> diff = new_slogan.fdiff(slogan)
    Diff(added=(Feature(FeatureLocation(ExactPosition(0), ExactPosition(9)), type='library'), Feature(FeatureLocation(ExactPosition(17), ExactPosition(18)), id='DNA')), removed=(Feature(FeatureLocation(ExactPosition(14), ExactPosition(17)), id='DNA'), Feature(FeatureLocation(ExactPosition(0), ExactPosition(5)), type='library'), Feature(FeatureLocation(ExactPosition(13), ExactPosition(16)), id='DNA')))
    >>> d.added
    (Feature(FeatureLocation(ExactPosition(0), ExactPosition(9)), type='library'),)
    >>> d.removed
    (Feature(FeatureLocation(ExactPosition(13), ExactPosition(16)), id='DNA'),
     Feature(FeatureLocation(ExactPosition(0), ExactPosition(5)), type='library'))

Internally, these values are stored in ``Component.features.added`` and ``Component.features.removed``.

.. note::

    Currently :meth:`Component.fdiff` is only implemented for components that directly inherit from one another.


Feature search
^^^^^^^^^^^^^^

- find method
- overlap method

Combining components
--------------------

Multiple components can be combined using :meth:`Component.combine`. This function will either create a `"source"`
feature annotation for each of the components that are being merged, or copy over all features from all components if
the keyword-only argument ``copy_features=True`` is set.

.. code-block:: python

    >>> a = Component('Co')
    >>> b = Component('Lib')
    >>> b.features.add(FeatureLocation(0, 3), id='lib')
    >>> c = Component.combine(a, b, copy_features=True)
    >>> c.seq
    Seq('CoLib', Alphabet())
    >>> c.features
    ComponentFeatureSet([Feature(FeatureLocation(ExactPosition(2), ExactPosition(5)), id='lib')])


Strain inheritance
------------------

In addition to DNA components, `colib` can track changes in haploid microbial organisms. :class:`HaploidOrganism`
can track added, changed, or deleted DNA components -- such as chromosomes or plasmids -- and aggregate features
contained in the strains.

- example
- find features in strain