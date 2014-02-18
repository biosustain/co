

class UnicellularOrganism(object):

    def __init__(self, display_id, parents=None):
        self.components = set()

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
