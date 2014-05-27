import weakref
import heapq
import logging
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys

from co.interval import IntervalTree, BaseInterval

class FeatureWrapper(BaseInterval):
    def __init__(self, feature):
        self.feature = feature

    @property
    def start(self):
        return self.feature.location.start

    @property
    def end(self):
        return self.feature.location.end


class FeatureSet(object):
    """
    An ordered collection of :class:`SeqFeature` objects.

    :param type feature_class: type of the features stored in the collection; defaults to :class:`SeqFeature` and must
        inherit from it.
    """

    def __init__(self, feature_class=None):
        if feature_class is None:
            feature_class = SeqFeature
        elif not issubclass(feature_class, SeqFeature):
            raise RuntimeError("FeatureSet expects a feature class that inherits from SeqFeature")

        self._features = IntervalTree()
        self._feature_class = feature_class

    def __or__(self, other):
        return self.difference(other)

    def __len__(self):
        return len(self._features)

    def __iter__(self):
        for f in self._features:
            if isinstance(f, FeatureWrapper):
                yield f.feature
            else:
                yield f

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, list(self))

    @staticmethod
    def _wrap_feature(feature):
        if isinstance(feature, BaseInterval):
            return feature
        else:
            return FeatureWrapper(feature)

    def copy(self):
        """
        :returns: a copy of this collection
        :rtype: :class:`FeatureSet`
        """
        fs = FeatureSet(feature_class=self._feature_class)
        fs._features = self._features.copy()
        return fs

    def add(self, *args, **kwargs):
        """
        Creates a feature object from the given ``args`` and ``kwargs`` and adds it to the collection.

        :rtype: :class:`SeqFeature`
        """
        feature = self._feature_class(*args, **kwargs)

        self._features.add(self._wrap_feature(feature))
        return feature

    def remove(self, feature):
        """
        Removes the given feature from the collection
        """
        self._features.remove(self._wrap_feature(feature))

    def find(self,
             between_start=None,
             between_end=None,
             type=None,
             id=None,
             strand=None,
             **qualifiers):
        """
        Iterate over all features matching the search parameters.

        - ``between_start`` and ``between_end`` can be used to restrict the search range.
        - ``type``, ``id``, and ``strand`` each restrict the search to features that match on these attributes
        - ``qualifiers`` is an arbitrary group of keyword arguments that will be matched to the qualifier keys of
          each feature. Each key must be present and have the same value as in the search parameters.

        """

        if between_start or between_end:
            it = self.overlap(between_start or 0, between_end or sys.maxsize)
        else:
            it = iter(self)

        attrs = [(k, v) for k, v in (('type', type), ('id', id), ('strand', strand)) if v is not None]

        for feature in it:
            if any(getattr(feature, key) != value for key, value in attrs):
                continue
            if any(feature.qualifiers.get(key) != value for key, value in qualifiers.items()):
                continue
            yield feature

    def overlap(self, start, end):
        """
        Returns an iterator over all features in the collection that overlap the given range.

        :param int start: overlap region start
        :param int end: overlap region end
        """
        if start > end:
            raise RuntimeError("start cannot be larger than end.")

        for f in self._features.find_overlapping(start, end):
            if isinstance(f, FeatureWrapper):
                yield f.feature
            else:
                yield f

    def difference(self, other):
        fs = self.copy()
        for f in other:
            fs._features.remove(f)
        return fs

    def union(self, other):
        fs = self.copy()
        for f in other:
            fs.add(f)
        return fs


class Feature(SeqFeature, BaseInterval):
    """
    :class:`Feature` derives from :class:`SeqFeature` and binds to a particular
    :class:`Component`. :class:`Feature` does not support the ``sub_features`` argument. All other
    :class:`SeqFeature` arguments are supported.
    """

    def __init__(self, component, *args, **kwargs):
        if 'sub_features' in kwargs:
            raise NotImplementedError('sub_features is not supported by ComponentSeqFeature and '
                                      'deprecated in SeqFeature. Use CompoundFeatureLocation instead.')

        self.source = kwargs.pop('source', None)
        self.component = component

        super(Feature, self).__init__(*args, **kwargs)

    def __hash__(self):
        return hash((self.location.start, self.location.end, self.type))

    def __eq__(self, other):
        return self.start == other.start and \
               self.end == other.end and \
               self.type == other.type and \
               self.strand == other.strand and \
               self.id == other.id and \
               self.qualifiers == other.qualifiers and \
               self.ref == other.ref and \
               self.ref_db == other.ref_db

    def _move_to_location(self, location, component=None):
        return Feature(location=location,
                       component=component or self.component,
                       type=self.type,
                       location_operator=self.location_operator,
                       id=self.id,
                       qualifiers=dict(self.qualifiers.items()))

        # TODO refs?

    def _move(self, start, end):
        # TODO ref
        new_location = FeatureLocation(start, end, strand=self.location.strand)
        return self._move_to_location(new_location)

    def _shift(self, offset, component=None):
        return self._move_to_location(self.location._shift(offset), component)

    @property
    def seq(self):
        """
        The sequence of the feature within the component as :class:`Seq` object.
        """
        seq = self.component.seq[self.start:self.end]
        if self.location.strand == REVERSE_STRAND:
            return seq.reverse_complement()
        return seq

    @property
    def start(self):
        """
        """
        return self.location.start

    @property
    def end(self):
        """
        """
        return self.location.end


class ComponentFeatureSet(FeatureSet):
    """
    An extended version of :class:`FeatureSet` that binds to a :class:`Component` and inherits from any
    :class:`FeatureSet` in the parent of a component.

    When iterating over this feature set, any inherited features are also returned.

    .. attribute:: removed_features

        Removed features are stored in :attr:`removed_features` if they are present in the parent, but not in
        :attr:`component`.

    .. attribute:: component

    """
    def __init__(self, component, removed_features=None, feature_class=None):
        super(ComponentFeatureSet, self).__init__(feature_class or Feature)
        self.component = component

        parent = component.parent
        self.parent_feature_set = parent.features if parent is not None else None
        self.removed_features = set(removed_features or ())

    def __iter__(self):
        if self.parent_feature_set:  # NOTE: this is where caching should kick in on any inherited implementation.
            keep_features = (f for f in self.parent_feature_set if f not in self.removed_features)
            translated_features = (self.component._translate_feature(f, self.component) for f in keep_features)
            return heapq.merge(self._features, translated_features)
        else:
            return super(ComponentFeatureSet, self).__iter__()

    def add(self, location, *args, **kwargs):
        if isinstance(location, self._feature_class):
            self._features.add(location)
            return location
        else:
            if self._feature_class == Feature:
                kwargs['component'] = self.component
            return super(ComponentFeatureSet, self).add(location=location, *args, **kwargs)

    @property
    def added(self):
        """
        Return the set of features added to this component, excluding any inherited features.
        """
        return self._features

    @property
    def removed(self):
        """
        Return all features present in the parent component's feature set but removed from this feature set.
        """
        return self.removed_features

    def remove(self, feature):
        if feature in self._features:
            self._features.remove(feature)
        elif self.parent_feature_set and feature in self.parent_feature_set:
            self.removed_features.add(feature)
        else:
            raise KeyError(feature)

    def overlap(self, start, end, include_inherited=True):
        """
        Returns an iterator over all features in the collection that overlap the given range.

        :param int start: overlap region start
        :param int end: overlap region end
        :param bool include_inherited: if ``True``, also yield all overlapping features in any ancestors of the
                                       component
        """
        intersect = set(super(ComponentFeatureSet, self).overlap(start, end))
        tt = self.component.tt()

        if self.parent_feature_set and include_inherited:
            if end >= tt.source_size:  # FIXME may be source_end
                end = tt.source_size - 1

            translated_start = tt.ge(start)
            translated_end = tt.le(end)

            logging.debug('find_overlapping({}, {}): '
                          'inherited between {} and {}'.format(start, end, translated_start, translated_end))

            if translated_start <= translated_end:
                intersect |= self.parent_feature_set.overlap(translated_start, translated_end) - self.removed_features
        return intersect


FORWARD_STRAND, REVERSE_STRAND = 1, -1


class Source(object):
    def __init__(self, component, broken=False):
        from co import Component
        if isinstance(component, Component):
            self._component_ref = weakref.ref(component)
            self._identifier = component.id
        else:
            self._identifier = component
        self.is_broken = broken

    @property
    def component(self):
        if self._component_ref:
            return self._component_ref()
        return None

    @property
    def identifier(self):
        return self._identifier

    def make_broken(self):
        if self.is_broken:
            return self
        return Source(self.component or self.identifier, broken=True)

    def __eq__(self, other):
        if isinstance(other, Source):
            return self.component == other.component and \
                   self.identifier == other.identifier and \
                   self.is_broken == other.is_broken
        return False

    def __repr__(self):
        return 'Source({}, broken={})'.format(repr(self.component or self.identifier), self.is_broken)