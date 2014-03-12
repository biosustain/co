from blinker import Signal


class HaploidOrganism(object):
    on_component_added = Signal()
    on_component_removed = Signal()

    def __init__(self, display_id, parents=None):
        self.components = set()

    @classmethod
    def from_file(cls, file, converter=None):
        pass

    def get_lineage(self):
        pass

    def is_locked(self):
        """
        An organism is locked for changes as soon as it has any children; an organism
        is locked for the creation of children while it is not registered in a storage.
        """
        pass

    def diff(self, other):
        """
        :param other:
        :type other: `HaploidOrganism`
        :returns: A list of added, changed, and removed component names.
        """
        pass
