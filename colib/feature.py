import weakref
import heapq
import logging
from Bio.SeqFeature import SeqFeature, FeatureLocation

from colib.interval import IntervalTree, BaseInterval

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

    def __init__(self, feature_class=None):
        if feature_class is None:
            feature_class = SeqFeature
        elif not issubclass(feature_class, (BaseInterval, SeqFeature)):
            raise RuntimeError("FeatureSet expects a feature class of type Interval or SeqFeature")

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

    def copy(self):
        fs = FeatureSet(feature_class=self._feature_class)
        fs._features = self._features.copy()
        return fs

    def add(self, *args, **kwargs):
        feature = self._feature_class(*args, **kwargs)

        if isinstance(feature, BaseInterval):
            self._features.add(feature)
        else:
            self._features.add(FeatureWrapper(feature))
        return feature

    def remove(self, feature):
        self._features.remove(feature)

    def find(self, between_start=None, between_end=None, **kwargs):
        raise NotImplementedError()

    def overlap(self, start, end):
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

    :class:`Feature` does not support the ``sub_features`` argument.

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
        seq = self.component.seq[self.start:self.end]
        if self.location.strand == REVERSE_STRAND:
            return seq.reverse_complement()
        return seq

    @property
    def start(self):
        return self.location.start

    @property
    def end(self):
        return self.location.end


class ComponentFeatureSet(FeatureSet):
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
        return self._features

    @property
    def removed(self):
        return self.removed_features

    def remove(self, feature):
        """
        """
        if feature in self._features:
            self._features.remove(feature)
        elif self.parent_feature_set and feature in self.parent_feature_set:
            self.removed_features.add(feature)
        else:
            raise KeyError(feature)

    def overlap(self, start, end, include_inherited=True):
        intersect = set(super(ComponentFeatureSet, self).overlap(start, end))
        tt = self.component.tt()

        if self.parent_feature_set and include_inherited:
            translated_start = tt.ge(start)
            translated_end = tt.le(end)

            logging.debug('find_overlapping({}, {}): '
                          'inherited between {} and {}'.format(start, end, translated_start, translated_end))

            if translated_start <= translated_end:
                intersect |= self.parent_feature_set.overlap(translated_start, translated_end) - self.removed_features
        return intersect


FORWARD_STRAND, REVERSE_STRAND = 1, -1
#
# class Feature(BaseInterval):
#     def __init__(self, component, position, size, strand=None, type=None, source=None, **attributes):
#         self._component = component
#         self._position = int(position)
#         self._size = int(size)
#         self._type = type
#
#         if 'source' in attributes:
#             source = attributes['source']
#             if not isinstance(source, Source):
#                 attributes['source'] = Source(source)
#
#         self.strand = strand
#         self._attributes = attributes
#
#     def translate_to(self, component, using_tt=None):
#         if component == self.component:
#             logging.warn('Translating {} from component {} to identical component'.format(self, self.component))
#             return self
#         if using_tt is None:
#             using_tt = component.tt(self.component)
#         print(self._position, using_tt[self.position])
#         return self.move(using_tt[self.position], size=self.size, component=component)
#
#     def move(self, position, size=None, component=None):
#         return self.__class__((component or self._component),
#                               position,
#                               size or self.size,
#                               strand=self.strand,
#                               type=self.type,
#                               **self._attributes)
#
#     @property
#     def type(self):
#         return self._type
#
#     @property
#     def qualifiers(self):
#         if 'qualifiers' in self._attributes:
#             return dict(self._attributes['qualifiers'])
#         return {}
#
#     @property
#     def component(self):
#         return self._component
#
#     @property
#     def position(self):
#         return self._position
#
#     @property
#     def size(self):
#         return self._size
#
#     @property
#     def component(self):
#         return self._component
#
#     @property
#     def start(self):
#         return self._position
#
#     @property
#     def end(self):
#         return self._position + self._size - 1
#
#     @property
#     def attributes(self):
#         return dict(self._attributes)
#
#     @property
#     def seq(self):
#         seq = self._component.seq[self.start:self.end + 1]
#         if self.strand == REVERSE_STRAND:
#             return seq.reverse_complement()
#         return seq
#
#     def is_anonymous(self):
#         return self._type is None
#
#     def is_unbound(self):
#         return self.component is None
#
#     def __eq__(self, other):
#         if not isinstance(other, Feature):
#             return False
#
#         return self.position == other.position and \
#                self.size == other.size and \
#                self.strand == other.strand and \
#                self.type == other.type and \
#                str(self.seq) == str(other.seq) and \
#                self.attributes == other.attributes
#
#     def __getattr__(self, item):
#         return self._attributes[item]
#
#     def __hash__(self):
#         return hash((self.start, self.end, self.strand, self.type, tuple(self.attributes)))
#
#     def __repr__(self):
#         args = '{}, {}, {}'.format(repr(self.component), self.position, self.size)
#         #return str(self.sequence)
#
#         if self.type:
#             args = '{}, type=\'{}\''.format(args, self.type)
#         if self.strand:
#             args = '{}, strand=\'{}\''.format(args, self.strand)
#         if self._attributes:
#             args = '{}, {}'.format(args, ', '.join('{}={}'.format(k, repr(v)) for k, v in self._attributes.items()))
#
#         return '{}({})'.format(self.__class__.__name__, args)
#


class Source(object):
    def __init__(self, component, broken=False):
        from colib import Component
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