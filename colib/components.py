# coding: utf-8
from collections import namedtuple
import functools
import heapq
import logging
from operator import attrgetter
from Bio.Seq import Seq
from blinker import Signal
import six
from colib.annotations import FeatureBase, FORWARD_STRAND
from colib.sequence import OverlapError, TranslationTable
from colib.utils import SortedCollection


def diff(original, other):
    if original.precursor == other:
        return original._mutations
    raise NotImplementedError()  # TODO if original.inherits(other)

class _FeatureList(object):

    def __init__(self, component, parent=None, feature_class=None):
        self._feature_class = feature_class or _Feature
        self._component = component
        self._features = SortedCollection(key=attrgetter('position'))
        self._features_from_end = SortedCollection(key=attrgetter('end'))
        self._removed_features = set()
        self._parent_feature_list = parent

    def intersect(self, start, end, include_inherited=True):
        if start > end:
            return set()

        # TODO this intersect lookup requires optimization.
        intersect = set(self._features_from_end.filter_ge(start)) & set(self._features.filter_le(end))

        if self._parent_feature_list and include_inherited:
            logging.debug('intersect({}, {}): inherited between {} and {}'.format(start, end, self._tt.ge(start), self._tt.le(end)))
            intersect |= self._parent_feature_list.intersect(self._tt.ge(start), self._tt.le(end)) - self._removed_features
        return intersect

    def find(self, *args, **kwargs):
        pass

    def add(self, position=None, size=None, **kwargs):
        if isinstance(position, FeatureBase):
            feature = position
            # feature = self._feature_class(self._component,
            #                    obj.position,
            #                    obj.size,
            #                    **obj.attributes)
        else:
            feature = self._feature_class(self._component, position, size, **kwargs)
        self._insert(feature)
        return feature

    def _insert(self, feature):
        self._features.insert(feature)
        self._features_from_end.insert(feature)

    def remove(self, feature):
        """
        """
        if feature in self._features:
            self._features.remove(feature)
            self._features_from_end.remove(feature)
        elif self._parent_feature_list:
            self._removed_features.add(feature)
        else:
            raise KeyError(feature)

    # def _translate_feature(self, feature):
    #     component = self._component
    #     return _Feature(component, component._mutations_tt[feature.position], *feature[2:])

    @property
    def _tt(self):
        return self._component._mutations_tt

    @property
    def added(self):
        return self._features

    @property
    def removed(self):
        return self._removed_features

    def __len__(self):
        if self._parent_feature_list:
            return len(self._parent_feature_list) - len(self._removed_features) + len(self._features)
        return len(self._features)

    def __iter__(self):
        if self._parent_feature_list:  # NOTE: this is where caching should kick in on any inherited implementation.
            keep_features = (f for f in self._parent_feature_list if f not in self._removed_features)
            translated_features = (f.translate(self._component, self._tt) for f in keep_features)

            for f in heapq.merge(self._features, translated_features):
                yield f
        else:
            for f in self._features:
                yield f


class Component(object):
    """
    If the underlying library supports it, mutated `Component` objects can be stored as a set of the features
    that have been added or removed from the parent component. It is recommended that the library implements a
    caching strategy for quickly looking up the precise list of features contained in a component.

    `Component` objects, once created, SHOULDN'T be edited. However they MAY be edited until they are referred
    to by another component -- either by direct mutation or through a `_Feature`. A strategy for deleting components
    without destroying descendant objects may be necessary.
    """
    on_feature_added = Signal()
    on_feature_removed = Signal()
    on_internal_mutation = Signal()

    def __init__(self, sequence='', parent=None, storage=None, meta=None, id=None, feature_class=None):
        if isinstance(sequence, six.string_types):
            sequence = Seq(sequence)

        # TODO support components without sequence.

        self._parent = parent
        self._sequence = sequence
        self._storage = storage
        self._mutations = ()
        self._mutations_tt = TranslationTable.from_mutations(self.sequence, self._mutations)
        self._id = id

        self.features = feature_list = _FeatureList(self,
                                                    feature_class=feature_class,
                                                    parent=self._parent.features if self._parent else None)
        self.meta = meta or {}



    @classmethod
    def from_file(cls, file, converter=None):
        pass

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
    def id(self):
        return self._id

    @property
    def length(self):
        return len(self._sequence)

    @property
    def parent(self):
        return self._parent

    @property
    def sequence(self):
        return self._sequence

    def clone(self):
        """
        Creates a copy of this :class:`Component`.

        The copy is linked to the same storage as the original, but is not saved automatically.
        """
        return Component(self._sequence,
                         parent=self,
                         storage=self._storage)

    def mutate(self, mutations, strict=True, clone=True):
        """
        Creates a copy of this :class:`Component` and applies all :param:`mutations`.

        Mutations are applied in order based on the original coordinate system. If `strict=True` and multiple
        mutations affect the same area, an `OverlapException` is raised. Otherwise, an `OverlapException` may
        still be raised if a mutation falls into a position that has previously been deleted.

        When a mutation affects a feature that has a linked component, that link will be flagged as `broken`.

        Memory
        ^^^^^^

        Although mutations are executed in order, the mutate function will attempt to resize features based on where
        their start or end based have been in the reference sequence if a based is lost and then restored in multiple
        mutations. For each deleted base, the algorithm tracks the next base to the left and to the right that was
        not deleted. In strict mode, since mutations are not allowed to overlap, this only affects directly adjacent
        mutations.

        ::

              ████      ████          ██        ██     ██        ██          ██       ████
            123  45   123  45    123  45   123  45    123  45   123  45    12345    123  45
            123xy45   123  -5    123xy45   123  -5    123xy45   123  -5    123-5    123xx45
            123xy-5   123xy-5    123xy-5   123xy-5    123xy-5   123xy-5      █░       ████
              ███░      ███░          ░█        ░█     ██        ██

        :param strict: If `True`, will fail if two mutations affect the same sequence.
        :param clone: `True` by default; if `False`, the component will *NOT* be cloned and the changes will be applied
            to the original version instead.
        :returns: A mutated copy of the component (or a reference to the original component if `clone=False`)
        :raises: OverlapException
        """
        component = self.clone() if clone else self
        features = _FeatureList(component, parent=self.features) if clone else self.features
        sequence = self._sequence.tomutable()

        tt = TranslationTable(self.length) if clone else self._mutations_tt
        changed_features = set()

        logging.debug('original features: {}'.format(list(self.features)))
        logging.debug('new features: {}'.format(list(features)))

        for mutation in mutations:
            # translate mutations, apply.
            # catch strict errors
            # flag features as edited/changed (copies are created as this happens), deleted.

            # FIXME handle addition to end of sequence.
            # FIXME requires proper handling of r/q_end to find out new size.

            logging.debug('---> mutation: "%s"', mutation)

            if strict:
                translated_start = tt[mutation.position]
            else:
                try:
                    translated_start = tt.le(mutation.position)
                    logging.debug('translated start: nonstrict=%s; strict=%s', translated_start, tt[mutation.position])
                except IndexError:
                    raise OverlapError()

            # TODO strict mode implementation that also fires on other overlaps.

            # TODO features.
            affected_features = features.intersect(mutation.position, mutation.position + mutation.size)


            logging.debug('affected features: {}'.format(affected_features))

            # TODO also search previously affected features and exclude these.
            # TODO e.g. features = _FeatureList(component, inherit=self.features)

            logging.debug('sequence: "%s"', sequence)
            logging.debug('alignment: \n%s', tt.alignment_str())

            if translated_start is None:
                raise OverlapError()

            logging.debug('new features: {}'.format(list(features)))
            logging.info('features in mutation: {}'.format(affected_features))

            for feature in affected_features | changed_features:
                # logging.debug(list(range(len(self.sequence))))
                # logging.debug(list(tt))
                logging.info('{} with sequence "{}" affected by {}.'.format(feature, feature.sequence, mutation))


                if feature.end < mutation.start or feature.start > mutation.end:
                    continue

                features.remove(feature)
                try:
                    changed_features.remove(feature)
                except KeyError:
                    pass

                # TODO find untranslated equivalents and replace when all is over.
                if mutation.start > feature.start:
                    if mutation.end < feature.end or mutation.size == 0:  # mutation properly contained in feature.
                        logging.debug('FMMF from {} to {}({})'.format(
                            feature.start, # tt.ge(feature.start),
                            feature.end, # tt.le(feature.end),
                            feature.end - mutation.size + mutation.new_size))

                        changed_features.add(feature.move(feature.start,
                                                          feature.size - mutation.size + mutation.new_size))
                    else:
                        logging.debug('FMFM from {} to {}'.format(
                            tt[feature.start],
                            tt[mutation.start - 1]))  # tt[mutation.start] ?

                        changed_features.add(feature.move(feature.start, mutation.start - feature.start))
                else:  # mutation.start >= feature.start

                    if mutation.end < feature.end:
                        logging.debug('MFMF from {} to {}'.format(tt[mutation.end + 1], tt[feature.end]))

                        changed_features.add(feature.move(mutation.end + 1, feature.end - mutation.end))
                    else:
                        pass # feature removed and not replaced.
                        logging.debug('MFFM feature removed')

            logging.debug('applying mutation: at %s("%s") translating from %s("%s") delete %s bases and insert "%s"',
                          translated_start,
                          sequence[translated_start],
                          mutation.start,
                          self.sequence[mutation.start],
                          mutation.size,
                          mutation.new_sequence)

            # apply mutation to sequence:
            if mutation.new_size != mutation.size:  # insertion, deletion or delins
                del sequence[translated_start:translated_start + mutation.size]
                tt.delete(translated_start, mutation.size)

                if mutation.new_size != 0:  # insertion or DELINS
                    sequence = sequence[:translated_start] + mutation.new_sequence + sequence[translated_start:]
                    tt.insert(translated_start, mutation.new_size)
            else:  # substitution:
                if mutation.size == 1:  # SNP
                    sequence[translated_start] = mutation.new_sequence
                else:
                    sequence = sequence[:translated_start]\
                        + mutation.new_sequence\
                        + sequence[translated_start + mutation.new_size:]

                tt.substitute(translated_start, mutation.size)

            logging.debug('new sequence: "%s"', sequence)
            logging.debug('changed features: {}'.format(changed_features))

        logging.debug('mutated sequence: %s', str(sequence))

        component._sequence = sequence.toseq()
        component._mutations = tuple(mutations)
        component._mutations_tt = tt
        component.features = features

        logging.debug('\n' + tt.alignment_str())

        for feature in changed_features:
            logging.debug('changed: %s', feature)
            translated_feature = feature.translate(component, tt)
            features._insert(translated_feature)

            logging.debug('translating %s: %s(%s) "%s" -> %s(%s) "%s"',
                          feature,
                          feature.position,
                          feature.size,
                          feature.sequence,
                          translated_feature.position,
                          translated_feature.size,
                          translated_feature.sequence)

        return component

    def get_lineage(self):
        component = self
        while self.parent:
            component = self.parent
            yield component

    def inherits(self, other):
        """
        Returns `True` if this object is a mutated version of `other`, `False` otherwise.
        :returns: Boolean
        """
        raise NotImplementedError()

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
        if other == self.parent:
            return self._mutations
        if self._storage and self._storage == other._storage:
            return self._storage.diff(self, other)
        return diff(self, other)


class _Feature(FeatureBase):

    def __init__(self, component, position, size, type=None, name=None, qualifiers=None, strand=FORWARD_STRAND, source=None, source_link_broken=False):
        super(_Feature, self).__init__(component, position, size, strand)
        self._type = type
        self._name = name
        self._qualifiers = qualifiers or {}
        self._source = source
        self._source_link_broken = source_link_broken

    component = property(fget=attrgetter('_component'))
    position = property(fget=attrgetter('_position'))
    size = property(fget=attrgetter('_size'))
    type = property(fget=attrgetter('_type'))
    name = property(fget=attrgetter('_name'))
    # strand = property(fget=attrgetter('_strand'))
    source = property(fget=attrgetter('_source'))
    source_link_broken = property(fget=attrgetter('_source_link_broken'))

    def translate(self, component, tt):
        # size = tt[self.end] - tt[self.start] + 1
        return self.move(tt[self.position], self.size, component=component)

    def move(self, position, size=None, component=None):
        return _Feature((component or self._component),
                        position,
                        size or self.size,
                        self._type,
                        self._name,
                        self._qualifiers,
                        self.strand,
                        self._source,
                        self._source_link_broken)

    @property
    def qualifiers(self):
        return self._qualifiers.copy()

    def is_anonymous(self):
        return self._type is None

    def __repr__(self):
        return '<Feature:{} from {} to {}>'.format(self.type or self.name, self.start, self.end)