from bisect import insort_left


class _StorageCollection(object):

    def __init__(self):
        self._item_list = []

    def add(self, value):
        insort_left(self._item_list, value)

    def bind(self): # TODO some sort of signaling to the library.
        pass

    def __getitem__(self, item):
        return self._item_list[item]

    def __iter__(self):
        return self._item_list.__iter__()


class NamedCollection(object):

    def __init__(self, name_getter=None):
        pass


class Storage(object):
    """


    A library takes care of loading and storing components. A component retrieved from the library may be lazily loaded.

    """
    components = _StorageCollection()
    organisms = _StorageCollection()

    def add_component(self, component):
        """
        The library implementation may also replace `component.features` with a collection implementation so that these
        do not have to be retrieved all at once when a component is looked up.
        """

        # Transparently converts component.features into a collection that may map to a database model.

    def add_organism(self, organism):
        pass

    def get_component(self, component):
        pass
