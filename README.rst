Co
==

**Co** is a Python library for altering annotated DNA sequences. It keeps track of components and lifts
over feature annotations when a component is "mutated" by applying a series of mutations. With ``co`` you can
build new consensus sequences for cloned organisms and follow changes to feature annotations within a lineage.

::

    >>> from co import Component
    >>> from co.mutations import *
    >>> hi_x = Component('Hello X!')
    >>> hi_x.seq
    Seq('Hello X!', Alphabet())
    >>> hi_world = hi_x.mutate([Mutation(6, 1, 'world')])
    >>> hi_world.seq
    Seq('Hello world!', Alphabet())

