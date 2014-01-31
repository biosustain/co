

class Annotation(object):
    """
    An annotation is different from a feature in that it is not bound to a particular `DNAComponent` and sequence.

    This means that annotations can be more easily applied to different components without being tightly coupled.

    Annotations are very simple objects almost like `Range` objects that are meant to be used by tools outside the
    library to attach additional information to components and being able to translate this information to mutated
    versions without having to deal with all the complexities of mutations.

    """
    def translate(self, from_component, to_component):
        pass

    def to_feature(self, component):
        pass


class AnnotationSet(object):

    def translate(self, from_component, to_component):
        pass