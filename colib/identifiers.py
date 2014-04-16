"""
TODO consider namedtuple.
"""

class UniqueIdentifier(object):
    def __init__(self, type, identifier): # TODO alt. names, databases
        self.type = type
        self.identifier = identifier

    def as_version(self, version_number):
        return Version(self.type, self.identifier, version_number)

    def __repr__(self):
        return '{}:{}'.format(self.type, self.identifier)


class Version(UniqueIdentifier):

    def __init__(self, type, identifier, version_number): # TODO alt. names, databases
        super(Version, self).__init__(type, identifier)
        self.version_number = version_number

    def __repr__(self):
        return '{}:{}.{}'.format(self.type, self.identifier, self.version_number)


class DatabaseXRef(UniqueIdentifier):
    pass
                     # /db_xref="GI:296146166"
                     # /db_xref="GeneID:853656"