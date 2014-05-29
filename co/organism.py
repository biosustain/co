import itertools
from co.difference import Diff


class FeatureView(object):
    """
    Iterates over all features in a set of components -- not necessarily in order -- and
    provides search access to these features.
    """

    def __init__(self, components):
        self.components = set(components)

    def find(self, **kwargs):
        """
        Searches all components and yields features matching the constraints.

        .. seealso:: :meth:`feature.FeatureSet.find`
        """
        return itertools.chain(*(c.find(**kwargs) for c in self.components))

    def __iter__(self):
        return itertools.chain(*(c.features for c in self.components))

    def __len__(self):
        return sum(len(c.features) for c in self.components)


class HaploidOrganism(object):
    """

    .. attribute:: id

        ID of the organism

    .. attribute:: parent

        Parent organism
    """

    def __init__(self, id, parent=None):
        self.id = id
        self.parent = parent
        self.components = {}

    def set(self, name, component):
        """

        :param str name: key for this component
            e.g. ``'genome'``, ``'chr1'``, or ``'pLASMID'``
        :param Component component:
        """
        self.components[name] = component

    def remove(self, name):
        """
        :param str name: key of the component to remove
        """
        del self.components[name]

    def get_lineage(self, inclusive=True):
        """
        Iterate over all ancestors of this organism

        :param bool inclusive: whether to include the organism itself in the lineage
        :returns: iterator over :class:`HaploidOrganism` objects
        """
        if inclusive:
            yield self

        for ancestor in self.parent.get_lineage():
            yield ancestor

    @property
    def features(self):
        """
        A read-only view of all the features present in all components of the organism.

        :returns: :class:`FeatureView`
        """
        return FeatureView(self.components.values())

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
            added=names - other_names,
            removed=other_names - names,
            changed=(name for name in names & other_names if self.components[name] != other.components[name]))

    def fdiff(self, other):
        """
        :returns: A :class:`Diff` object with added and removed features.
        """
        pass