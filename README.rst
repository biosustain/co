Co
==

**Co** is a Python library for altering annotated DNA sequences. It keeps track of components and lifts
over feature annotations when a component is "mutated" by applying a series of mutations. With ``co`` you can
build new consensus sequences for cloned organisms and follow changes to feature annotations within a lineage.

For more information, check out the `Documentation <http://co.readthedocs.org/en/latest/>`_.

Hello Co!
---------

::

    >>> from co import Component
    >>> from co.mutation import *
    >>> hello = Component('Hello X!')
    >>> hello.seq
    Seq('Hello X!', Alphabet())
    >>> hello_world = hello.mutate([Mutation(6, 1, 'world')])
    >>> hello_world.seq
    Seq('Hello world!', Alphabet())



Working with Feature Annotations
--------------------------------

Components are modeled after BioPython's ``SeqRecord`` -- they have both a sequence, and features:

.. code-block:: python

    >>> from Bio.SeqFeature import *
    >>> slogan = Component('CoPy is for DNA components', features=[
    ...                 SeqFeature(FeatureLocation(0, 4), type='name'),
    ...                 SeqFeature(FeatureLocation(12, 15), id='DNA')])
    >>>
    >>> # features are bound to components -- and you can always access their DNA sequence
    ...
    >>> slogan.features.add(FeatureLocation(16, 26)).seq
    Seq('components', Alphabet())
    >>> [f.seq for f in slogan.features]
    [Seq('CoPy', Alphabet()), Seq('DNA', Alphabet()), Seq('components', Alphabet())]
    >>>
    >>> # New Components are made through series of mutations
    ... # You not only get the new sequence but a mutated component: Features are translated to the
    ... # new sequence as well.
    ...
    >>> new_slogan = slogan.mutate([DEL(2, 2), DEL(12, 4)])
    >>> new_slogan.seq
    Seq('Co is for components', Alphabet())
    >>> new_slogan.features
    ComponentFeatureSet([Feature(FeatureLocation(ExactPosition(0), ExactPosition(2)), type='name'),
                         Feature(FeatureLocation(ExactPosition(10), ExactPosition(20)))])
    >>> [f.seq for f in new_slogan.features]
    [Seq('Co', Alphabet()), Seq('components', Alphabet())]
    >>> list(new_slogan.features.find(type='name'))  # features can be filtered by type, id, strand, position, and qualifiers
    [Feature(FeatureLocation(ExactPosition(0), ExactPosition(2)), type='name')]
    >>>
    >>> # Using Component.fdiff you can get a summary of what features where affected by mutation. (Unchanged features
    ... # that have a new coordinate -- e.g. the 'components' feature in this example -- are not included).
    ...
    >>> slogan.fdiff(new_slogan)
    Diff(added=(Feature(FeatureLocation(ExactPosition(12), ExactPosition(15)), id='DNA'),
                Feature(FeatureLocation(ExactPosition(0), ExactPosition(4)), type='name')),
         removed=(Feature(FeatureLocation(ExactPosition(0), ExactPosition(2)), type='name'),))


Authors
=======

`Lars Schöning <https://github.com/lyschoening>`_ has created Co. Contributions are very welcome.
Contact the main author for bigger changes.

