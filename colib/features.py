import heapq
import logging
from colib.interval import IntervalTree, IntervalBase


class FeatureSet(object):

    def __init__(self, feature_class=None):
        if feature_class is None:
            feature_class = Feature
        elif not issubclass(feature_class, IntervalBase):
            raise RuntimeError("FeatureSet expects a feature class of type Interval")

        self._features = IntervalTree()
        self._feature_class = feature_class

    def __or__(self, other):
        return self.difference(other)

    def __len__(self):
        return len(self._features)

    def __iter__(self):
        for f in self._features:
            yield f

    def copy(self):
        fs = FeatureSet(feature_class=self._feature_class)
        fs._features = self._features.copy()
        return fs

    def add(self, *args, **kwargs):
        feature = self._feature_class(*args, **kwargs)
        self._features.add(feature)
        return feature

    def remove(self, feature):
        self._features.remove(feature)

    def find(self, between_start=None, between_end=None, **kwargs):
        raise NotImplementedError()

    def overlap(self, start, end):
        if start > end:
            raise RuntimeError("start cannot be larger than end.")

        return self._features.find_overlapping(start, end)

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


class ComponentFeatureSet(FeatureSet):

    def __init__(self, component, **kwargs):
        super(ComponentFeatureSet, self).__init__(**kwargs)
        self.component = component

        parent = component.parent
        self.parent_feature_set = parent.features if parent is not None else None
        self.removed_features = set()

    def __iter__(self):
        if self.parent_feature_set:  # NOTE: this is where caching should kick in on any inherited implementation.
            keep_features = (f for f in self.parent_feature_set if f not in self.removed_features)
            translated_features = (f.translate_to(self.component) for f in keep_features)

            for f in heapq.merge(self._features, translated_features):
                yield f
        else:
            for f in self._features:
                yield f

    def add(self, position, *args, **kwargs):
        if isinstance(position, self._feature_class):
            self._features.add(position)
            return position
        else:
            return super(ComponentFeatureSet, self).add(self.component, position, *args, **kwargs)

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
        elif self.parent_feature_set:
            self.removed_features.add(feature)
        else:
            raise KeyError(feature)

    def overlap(self, start, end, include_inherited=True):
        intersect = set(super(ComponentFeatureSet, self).overlap(start, end))
        tt = self.component.tt()

        if self.parent_feature_set and include_inherited:
            logging.debug('find_overlapping({}, {}): inherited between {} and {}'.format(start, end, tt.ge(start), tt.le(end)))
            intersect |= self.parent_feature_set.overlap(tt.ge(start), tt.le(end)) - self.removed_features
        return intersect


FORWARD_STRAND, REVERSE_STRAND = 1, -1


class Feature(IntervalBase):

    def __init__(self, component, position, size, strand=None, type=None, **attributes):
        self._component = component
        self._position = int(position)
        self._size = int(size)
        self._type = type
        self.strand = strand
        self._attributes = attributes

    def translate_to(self, component, using_tt=None):
        if component == self.component:
            logging.warn('Translating {} from component {} to identical component'.format(self, self.component))
            return self
        if using_tt is None:
            using_tt = component.tt(self.component)
        print(self._position, using_tt[self.position])
        return self.move(using_tt[self.position], size=self.size, component=component)

    def move(self, position, size=None, component=None):
        return Feature((component or self._component),
                       position,
                       size or self.size,
                       strand=self.strand,
                       type=self.type,
                       **self._attributes)

    @property
    def type(self):
        return self._type

    @property
    def qualifiers(self):
        if 'qualifiers' in self._attributes:
            return dict(self._attributes['qualifiers'])
        return {}

    @property
    def component(self):
        return self._component

    @property
    def position(self):
        return self._position

    @property
    def size(self):
        return self._size

    @property
    def component(self):
        return self._component

    @property
    def start(self):
        return self._position

    @property
    def end(self):
        return self._position + self._size - 1

    @property
    def sequence(self):
        sequence = self._component[self.start:self.end + 1]
        if self.strand == REVERSE_STRAND:
            return sequence.reverse_complement()
        return sequence

    def is_anonymous(self):
        return self._type is None

    def is_unbound(self):
        return self.component is None

    def __eq__(self, other):
        if not isinstance(other, Feature):
            return False
        return self.start == other.start and \
               self.size == other.size and \
               self.type == other.type and \
               self.strand == other.strand

    def __getattr__(self, item):
        return self._attributes[item]

    def __hash__(self):
        return hash((self.start, self.end, self.type))

    def __repr__(self):
        return '<Feature:{} from {} to {}>'.format(self.type or self.name, self.start, self.end)


class Annotation(object):
    """
    An annotation is different from a feature in that it is not bound to a particular `Component` and sequence.

    This means that annotations can be more easily applied to different components without being tightly coupled.

    Annotations are very simple objects almost like `Range` objects that are meant to be used by tools outside the
    library to attach additional information to components and being able to translate_to this information to mutated
    versions without having to deal with all the complexities of mutations.

    """
    def __init__(self, position, size, **attributes):
        self.position = position
        self.size = size
        self._attributes = attributes

    def __getattr__(self, item):
        return self._attributes[item]

    def translate(self, from_component, to_component):
        raise NotImplementedError()

