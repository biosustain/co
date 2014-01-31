from Bio.Seq import Seq
import six
from colib.mutations import TranslationTable
from colib.sequence import Sequence


def diff(original, other):
    pass


class Component(Sequence):
    """


    If the underlying library supports it, mutated `Component` objects can be stored as a set of the features
    that have been added or removed from the parent component. It is recommended that the library implements a
    caching strategy for quickly looking up the precise list of features contained in a component.

    `Component` objects, once created, SHOULDN'T be edited. However they MAY be edited until they are referred
    to by another component -- either by direct mutation or through a `_Feature`. A strategy for deleting components
    without destroying descendant objects may be necessary.
    """

    def __init__(self, sequence='', storage=None):
        if isinstance(sequence, six.string_types):
            sequence = Seq(sequence)

        # TODO support components without sequence.

        self._sequence = sequence
        self._storage = storage
        self.features = set()

    @property
    def length(self):
        return len(self._sequence)

    @property
    def sequence(self):
        return self._sequence

    def clone(self):
        """
        Creates a copy of this :class:`Component`.

        The copy is linked to the same storage as the original, but is not saved automatically.
        """
        pass

    def mutate(self, mutations, strict=True, clone=True):
        """
        Creates a copy of this :class:`Component` and applies all :param:`mutations`.

        Mutations are applied relative to the original coordinate systems.

        :param strict: If `True`, will fail if two mutations affect the same sequence.
        :param clone: `True` by default; if `False`, the component will *NOT* be cloned and the changes will be applied
            to the original version instead.
        """
        # """
        #
        # Strategies apply to how mutations to features are dealt with.
        #
        #     When an insertion or deletion overlaps the beginning or end of a feature:
        #
        #     - Keep the size of the feature the same.
        #     - Attempt to keep the size of the feature identical if no more than one end of the feature is
        #       affected and the adjacent feature is anonymous (i.e. carries no information).
        #
        #     When a mutation is applied to a feature that has an associated part:
        #
        #     - Mutate the part and create a copy.
        #     - Mutate the feature to an inexact mapping.
        #
        # Mutations are applied in order based on the original coordinate system. If `strict=True` and multiple
        # mutations affect the same area, an `OverlapError` is raised.
        #
        # Anonymous features are merged by default.
        #
        # :param mutations: The list of mutations to apply (in order)
        # :type mutations:
        # :returns: A mutated copy of the component and a list of any new components created.
        # :raises: OverlapError
        #
        # """

        component = self.clone() if clone else self
        sequence = self._sequence.tomutable()
        tt = TranslationTable(self.length)

        for mutation in mutations:
            # translate mutations, apply.
            # catch strict errors
            # flag features as edited/changed (copies are created as this happens), deleted.

            # TODO check if mutation is within original table.

            # deletion if new_sequence == None
            # deletion, followed by insertion if size > 1
            if mutation.new_sequence is None or mutation.size != 1:
                tt.delete(tt[mutation.start], mutation.size)

            if len(mutation.new_sequence) > 1:  # insertion
                # FIXME raises TypeError if mutation.start is out of bounds
                # TODO unit tests for edge cases
                tt.insert(tt[mutation.start] + 1, len(mutation.new_sequence) - 1)
            else:  # substitution when size == 1 and len(new_sequence) == 1
                tt.substitute(tt[mutation.start], 1)


            # TODO sequence[] =
            # TODO features ..

    def inherits(self, other):
        """
        Returns `True` if this object is a mutated version of `other`, `False` otherwise.
        :returns: Boolean
        """
        raise NotImplementedError

    def diff(self, other):
        """
        Returns the list of mutations that represent the difference between this component and another.

        .. note::

            This may be aided by the storage library. The storage library (if attached) should be called first to see if
            there is some sort of optimized lookup. The library itself will either implement an optimized version or
            fall back to the diff solution used here.
        """
        assert isinstance(other, Component)
        if self._storage and self._storage == other._storage:
            return self._storage.diff(self, other)
        return diff(self, other)


class _Feature(Sequence):
    def __init__(self, component, start, end, type, properties):
        pass

    def is_anonymous(self):
        pass

    def identical_to(self, other):
        """
        Two features are identical if their sequence (and therefore length) is identical and their annotation
        is the same.
        """
        raise NotImplementedError