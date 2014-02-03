from collections import namedtuple
from Bio.Seq import Seq
import six
from colib.mutations import TranslationTable
from colib.sequence import Sequence


def diff(original, other):
    pass


class OverlapException(Exception):
    pass


class _FeatureList(object):

    def __init__(self, component, inherit=None):
        self._component = component
        self._features = set()
        self._removed_features = set()
        self._inherits = inherit

    def intersect(self, start, end, include_inherited=True):
        pass

    def find(self, *args, **kwargs):
        pass

    def add(self, position, size, type=None, name=None, attributes=None, link=None):
        feature = _Feature(self._component, position, size, type, name, attributes, link)
        self._features.add(feature)
        return feature

    def remove(self, feature):
        """
        """
        self._removed_features.add(feature)

    def __iter__(self):
        # sorted iteration over this list and the inherited list, automatically translating feature
        # positions from any inherited tables, excluding any in removed_features list.
        # TODO cache support for feature ids & positions.
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

    def __init__(self, sequence='', ancestor=None, storage=None):
        if isinstance(sequence, six.string_types):
            sequence = Seq(sequence)

        # TODO support components without sequence.

        self._ancestor = ancestor
        self._sequence = sequence
        self._storage = storage
        self.features = _FeatureList(self, inherit=self._ancestor.features if self._ancestor else None)

    @classmethod
    def from_components(cls, components, copy_features=False):
        """
        Joins multiple components together, creating a feature annotation for each.

        .. warning::

            Optionally, features from the old components can be copied into the new component. This is strongly
            discouraged when working with large components and should only be done when assembling a component from
            temporary components that will later be discarded. A custom inspection can recursively enumerate features
            from linked components.

        :param copy_features: When `copy_features` is set, no features with links to the original components will be
            created.
        """
        if copy_features:
            raise NotImplementedError("Copying features from old components into new ones is not currently supported.")

        # TODO

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
        return Component(self._sequence,
                         ancestor=self,
                         storage=self._storage)


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
        # mutations affect the same area, an `OverlapException` is raised.
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

        resolved_mutations = map(lambda m: m.resolve(self.sequence), mutations)

        for mutation in resolved_mutations:
            # translate mutations, apply.
            # catch strict errors
            # flag features as edited/changed (copies are created as this happens), deleted.

            insert_size = len(mutation.new_sequence)
            translated_start = tt[mutation.start]

            # TODO strict mode implementation that also fires on other overlaps.

            # TODO features.

            if translated_start is None:
                raise OverlapException()

            if insert_size != mutation.size:
                del sequence[translated_start:translated_start + mutation.size]
                tt.delete(translated_start, mutation.size)

                if insert_size != 0:  # insertion or DELINS
                    sequence = sequence[:translated_start] + mutation.new_sequence + sequence[translated_start:]
                    tt.insert(translated_start, insert_size)

            else:  # substitution
                if mutation.size == 1:  # SNP
                    sequence[translated_start] = mutation.new_sequence
                else:
                    sequence = sequence[:translated_start]\
                        + mutation.new_sequence\
                        + sequence[translated_start + insert_size:]

                tt.substitute(translated_start, mutation.size)

        # TODO features

        component._sequence = sequence.toseq()
        component._ancestor_mutations = tuple(resolved_mutations)
        return component

    def inherits(self, other):
        """
        Returns `True` if this object is a mutated version of `other`, `False` otherwise.
        :returns: Boolean
        """
        raise NotImplementedError

    # noinspection PyProtectedMember
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


class _Feature(namedtuple('Feature', ['component', 'position', 'size', 'type', 'name', 'attributes', 'link'])):
    __slots__ = ()

    @property
    def is_anonymous(self):
        return self.type is None

    def __gt__(self, other):
        return self.position > other.position
