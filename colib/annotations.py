FORWARD_STRAND, REVERSE_STRAND = 1, -1


class FeatureBase(object):

    def __init__(self, component, position, size, strand=None):
        self._component = component
        self._position = int(position)
        self._size = int(size)
        self.strand = strand

    def translate(self, component, tt):
        # size = tt[self.end] - tt[self.start] + 1
        return self.move(tt[self.position], self.size, component=component)

    def move(self, position, size, component=None):
        return FeatureBase(component or self._component, position, size,
                           strand=self.strand)

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
        sequence = self._component.sequence[self.start:self.end + 1]
        if self.strand == REVERSE_STRAND:
            return sequence.reverse_complement()
        return sequence

    def __lt__(self, other):
        return self._position < other._position

    def __eq__(self, other):
        return self._position == other._position

    def __hash__(self):
        return self._position


class Annotation(object):
    """
    An annotation is different from a feature in that it is not bound to a particular `Component` and sequence.

    This means that annotations can be more easily applied to different components without being tightly coupled.

    Annotations are very simple objects almost like `Range` objects that are meant to be used by tools outside the
    library to attach additional information to components and being able to translate this information to mutated
    versions without having to deal with all the complexities of mutations.

    """
    def __init__(self, position, size, **attributes):
        self.position = position
        self.size = size
        self._attributes = attributes

    @property
    def attributes(self):
        return self.attributes.copy()

    def translate(self, from_component, to_component):
        pass