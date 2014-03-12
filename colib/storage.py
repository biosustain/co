from bisect import insort_left
from colib.components import Component
from colib.organisms import HaploidOrganism


class _MemoryStorageCollection(object):

    def __init__(self, storage):
        self.storage = storage
        self._item_list = []

    def add(self, value):
        insort_left(self._item_list, value)

    def bind(self): # TODO some sort of signaling to the library.
        pass

    def __getitem__(self, item):
        return self._item_list[item]

    def __iter__(self):
        return self._item_list.__iter__()


class ComponentStore(_MemoryStorageCollection):

    def add(self, value):
        pass


class FeatureStore(_MemoryStorageCollection):

    def filter(self, name=None, type=None, id=None):
        pass


class NamedCollection(object):

    def __init__(self, name_getter=None):
        pass


class Storage(object):
    """


    A library takes care of loading and storing components. A component retrieved from the library may be lazily loaded.

    """

    def __init__(self):
        self.components = _MemoryStorageCollection(self)
        self.organisms = _MemoryStorageCollection(self)

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


class _StorageResource(object):
    def __init__(self, storage, type):
        self._storage = storage
        self._resource_type = type

    def add(self, item):
        self._storage._add_item(item, self._resource_type)

    def remove(self, item):
        self._storage._remove_item(item, self._resource_type)

    def get(self, item_id):
        pass

    def find(self, **query):
        self._storage._add_item(self._resource_type, **query)


class BaseStorage(object):
    def __init__(self):
        self.components = _StorageResource(self, Component)
        self.organisms = _StorageResource(self, HaploidOrganism)


class RESTStorage(BaseStorage):
    def __init__(self, base_uri, client_id=None, token=None):
        super(RESTStorage, self).__init__()
        # self._base_uri = base_uri
        # self._session = OAuth2Session()

    def post(self, uri, data, **kwargs):
        return self._session.post(self._base_uri + uri, data=data, **kwargs)

    def get(self, uri, **kwargs):
        return self._session.get(self._base_uri + uri, **kwargs)

    def add_component(self, component):
        self.post('/components', data=None)

    def get_component(self, component_id):
        data = self.get('/components/{}'.format(component_id))



    # def remove_item(self, item, resource_type=None):
    #     pass
    #
    # def find_item(self, resource_type, **query):
    #     pass


class MemoryStorage(object):
    def __init__(self):
        pass
