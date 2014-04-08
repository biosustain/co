"""
TODO consider namedtuple.
"""

class UniqueIdentifier(object):
    def __init__(self, type, identifier): # TODO alt. names, databases
        self.type = type
        self.identifier = identifier


class Version(object):
    pass


class DatabaseXRef(UniqueIdentifier):
    pass
                     # /db_xref="GI:296146166"
                     # /db_xref="GeneID:853656"