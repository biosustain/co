
==========
Quickstart
==========

.. module:: co


Simple mutation
---------------


To illustrate what :mod:`co` is designed for, let's begin with a hello world example:

.. code-block:: python

    >>> from co import Component
    >>> from co.mutation import *
    >>> hello = Component('Hello X!')
    >>> hello.seq
    Seq('Hello X!', Alphabet())
    >>> hello_world = hello.mutate([Mutation(6, 1, 'world')])
    >>> hello_world.seq
    Seq('Hello world!', Alphabet())


.. module:: co
.. features

Features & feature inheritance
------------------------------

:attr:`Component.feature` stores feature annotations in :class:`Feature` format. These feature annotations are attached
to the component they are defined in allowing, among other things, easy access to a feature sequence. Features are
attached to a component like so:

.. code-block:: python

    >>> from Bio.SeqFeature import *
    >>>
    >>> slogan = Component('CoPy is for DNA components', features=[
    ...                 SeqFeature(FeatureLocation(0, 4), type='name'),
    ...                 SeqFeature(FeatureLocation(12, 15), id='DNA')])
    >>>
    >>> slogan.features.add(FeatureLocation(16, 26)).seq
    Seq('components', Alphabet())
    >>> [f.seq for f in slogan.features]
    [Seq('CoPy', Alphabet()), Seq('DNA', Alphabet()), Seq('components', Alphabet())]

When a component is mutated, :mod:`co` automatically translates the feature annotations from the parent to
the new coordinate system:

.. code-block:: python

    >>> new_slogan = slogan.mutate([DEL(2, 2), DEL(12, 4)])
    >>> new_slogan.seq
    Seq('Co is for components', Alphabet())
    >>> new_slogan.features
    ComponentFeatureSet([Feature(FeatureLocation(ExactPosition(0), ExactPosition(2)), type='name'),
                         Feature(FeatureLocation(ExactPosition(10), ExactPosition(20)))])
    >>>
    >>> [f.seq for f in new_slogan.features]
    [Seq('Co', Alphabet()), Seq('components', Alphabet())]


When a region is affected by a mutation, any features contained in that region are deleted. Features that overlap
the mutated region are trimmed. Features containing mutations are marked as changed. Features that are not affected by any
mutation are left as they were---their starting coordinates are rewritten on the fly to map to the coordinate system
of any inheriting component.

In the sample above, the `type` ``'name'`` feature is resized by the mutation. The sequence of the
`id` ``'DNA'`` feature is deleted in its entirety and so the feature is deleted too. The feature spanning ``'components'``
has not changed at all---but the mutations do affect its coordinates and so they are lifted over when the feature
is accessed from within the mutated component.

.. code-block:: python

    >>> new_slogan.features.removed
    {Feature(FeatureLocation(ExactPosition(0), ExactPosition(9)), type='name'),
     Feature(FeatureLocation(ExactPosition(17), ExactPosition(20)), id='DNA')}
    >>> list(new_slogan.features.added)
    [Feature(FeatureLocation(ExactPosition(0), ExactPosition(5)), type='name')]


Feature diffs
^^^^^^^^^^^^^

:meth:`Component.fdiff` is designed for comparing the features contained in any two components:

.. code-block:: python

    >>> diff = new_slogan.fdiff(slogan)
    Diff(added=(Feature(FeatureLocation(ExactPosition(0), ExactPosition(9)), type='library'), Feature(FeatureLocation(ExactPosition(17), ExactPosition(18)), id='DNA')), removed=(Feature(FeatureLocation(ExactPosition(14), ExactPosition(17)), id='DNA'), Feature(FeatureLocation(ExactPosition(0), ExactPosition(5)), type='library'), Feature(FeatureLocation(ExactPosition(13), ExactPosition(16)), id='DNA')))
    >>> d.added
    (Feature(FeatureLocation(ExactPosition(0), ExactPosition(9)), type='library'),)
    >>> d.removed
    (Feature(FeatureLocation(ExactPosition(13), ExactPosition(16)), id='DNA'),
     Feature(FeatureLocation(ExactPosition(0), ExactPosition(5)), type='library'))



.. note::

    :meth:`Component.fdiff` is currently only implemented for components that directly inherit from one another.
    Internally, these values are looked up from ``Component.features.added`` and ``Component.features.removed``
    as shown earlier. Eventually this will work with any two components regardless of ancestry.

Feature search
^^^^^^^^^^^^^^

Features can be filtered using :meth:`FeatureSet.find`. This search function supports filtering by region, type, id,
strand as well as any qualifier. Multiple search parameters are interpreted as logical "AND"---i.e. all of them have
to match.

.. code-block:: python

    >>> from co import *
    >>> from Bio.SeqFeature import *
    >>>
    >>> letters = Component('AABBDDEE', features=[
    ...             SeqFeature(FeatureLocation(0, 1), type='vowel'),
    ...             SeqFeature(FeatureLocation(2, 5), type='consonant'),
    ...             SeqFeature(FeatureLocation(5, 6), type='vowel', qualifiers={'gene': 'abcD'})])
    >>>
    >>> list(letters.features.find(type='vowel'))
    [Feature(FeatureLocation(ExactPosition(0), ExactPosition(1)), type='vowel'), Feature(FeatureLocation(ExactPosition(5), ExactPosition(6)), type='vowel')]
    >>> list(letters.features.find(between_start=3))
    [Feature(FeatureLocation(ExactPosition(5), ExactPosition(6)), type='vowel'), Feature(FeatureLocation(ExactPosition(2), ExactPosition(5)), type='consonant')]
    >>>
    >>> from co.mutation import *
    >>> letters = letters.mutate([INS(4, 'CC')])
    >>> letters.seq
    Seq('AABBCCDDEE', Alphabet())
    >>> list(letters.features.find(type='consonant'))
    [Feature(FeatureLocation(ExactPosition(2), ExactPosition(7)), type='consonant')]
    >>> list(letters.features.find(type='vowel'))
    [Feature(FeatureLocation(ExactPosition(0), ExactPosition(1)), type='vowel'), Feature(FeatureLocation(ExactPosition(7), ExactPosition(8)), type='vowel')]
    >>> list(letters.features.find(type='consonant', between_end=1))
    []
    >>> list(letters.features.find(gene='abcD'))
    [Feature(FeatureLocation(ExactPosition(7), ExactPosition(8)), type='vowel')]


Optimization behind the scenes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Feature annotations that are inherited from another component are not copied over
in memory --- instead they are looked up on the fly. Only added and removed features are stored. A feature is
considered changed when its sequence is affected in any way. When a feature is changed, the old feature is removed and
the new feature is added.

- On-the-fly coordinate translation of unchanged features is done using :class:`translation.TranslationTable`---inspired
  by the UCSC Chain Format.
- Feature locations are indexed using :class:`interval.IntervalTree`, currently implemented as a BST.


Combining components
--------------------

Multiple components can be combined using :meth:`Component.combine`. This function will either create a `"source"`
feature annotation for each of the components that are being merged, or copy over all features from all components if
``copy_features=True``.

.. code-block:: python

    >>> a = Component('Co')
    >>> b = Component('Py')
    >>> b.features.add(FeatureLocation(0, 3), id='python')
    >>> c = Component.combine(a, b, copy_features=True)
    >>> c.seq
    Seq('CoPy', Alphabet())
    >>> c.features
    ComponentFeatureSet([Feature(FeatureLocation(ExactPosition(2), ExactPosition(5)), id='python')])


Strain inheritance
------------------

In addition to DNA components, `co` can track changes in haploid microbial organisms. :class:`HaploidOrganism`
can track added, changed, or deleted DNA components---such as chromosomes or plasmids---and aggregate features
contained in the strains.

Strain components
^^^^^^^^^^^^^^^^^

:meth:`HaploidOrganism.diff` tracks how components have changed across strains:

    >>> from co.organism import *
    >>> from co import *
    >>>
    >>> genome = Component('A')
    >>> alpha = HaploidOrganism('alpha')
    >>> alpha.set('genome', genome)
    >>>
    >>> beta = HaploidOrganism('beta', parent=alpha)
    >>> beta.set('genome', genome.mutate([Mutation(0, 1, 'B')]))
    >>> beta.set('plasmid', Component('AGCT'))
    >>> beta.diff(alpha)
    Diff(added=(), removed=('plasmid',), changed=('genome',))
    >>> ~beta.diff(alpha)
    Diff(added=('plasmid',), removed=(), changed=('genome',))


Strain features
^^^^^^^^^^^^^^^

:attr:`HaploidOrganism.features` returns a :class:`organism.FeatureView` which is a searchable and iterable
view of all features in all components of a strain.
