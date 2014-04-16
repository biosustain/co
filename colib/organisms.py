
from colib.difference import Diff





class HaploidOrganism(object):
    """

    .. attribute::

    """
    def __init__(self, display_id, parent=None):
        self.display_id = display_id
        self.parent = parent
        self.components = {}
        self.component_types = {}

    def add(self, name, component, type_=None):
        self.components[name] = component
        self.component_types = type_

    def remove(self, name):
        del self.components[name]
        del self.component_types[name]

    def get_lineage(self, inclusive=True):
        if inclusive:
            yield self

        for ancestor in self.parent.get_lineage():
            yield ancestor

    def is_locked(self):
        """
        An organism is locked for changes as soon as it has any children; an organism
        is locked for the creation of children while it is not registered in a context.
        """
        raise NotImplementedError()

    @property
    def features(self):
        """
        A read-only view of all the features present in all components of the organism.
        """
        raise NotImplementedError()

    def diff(self, other):
        """
        :param other:
        :type other: `HaploidOrganism`
        :returns: A :class:`Diff` object with added, changed, and removed component names.
        """
        assert isinstance(other, HaploidOrganism)

        names = set(self.components.keys())
        other_names = set(other.components.keys())

        return Diff(
            added=other_names - names,
            removed=names - other_names,
            changed=(name for name in names & other_names if self.components[name] != other.components[name]))

    def fdiff(self, other):
        """
        :returns: A :class:`Diff` object with added and removed features.
        """
        pass