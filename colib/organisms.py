from bisect import bisect_left, bisect_right
from colib.position import Position


class UnicellularOrganism(object):

    def __init__(self, display_id, parents=None):
        self.contigs = set()

    @classmethod
    def from_file(cls, file, converter=None):
        pass

    def get_lineage(self):
        pass

    def diff(self, other):
        """
        :param other:
        :type other: `UnicellularOrganism`
        :returns: A list of added, changed, and removed `Contig` objects.
        """
        pass


class Contig(object):
    def __init__(self, name, component):
        pass


class FeatureCollection(object):
    """

    To be efficiently looked up, the following rules apply to features:

    - features are in order
    - features may overlap
    - features may not contain other features (i.e. overlap them by more bases than they are long)
    - when a mutation deletes all parts of a feature that do not overlap with the previous feature, that feature must
      be removed in its entirety.

    """

    def __init__(self):
        self._features = []

    def create(self, start, end, sequence, type=None, meta=None):
        """
        :raises: OverlapException
        """

    def find(self, start, end):
        """
        """
        start_idx = bisect_left(self._features, start)
        end_idx = bisect_right(self._features, end)

        return [self._features[i] for i in range(start_idx, end_idx)]

    def __iter__(self):
        pass