"""
TODO consider namedtuple.
"""


class UniqueIdentifier(object):
    def __init__(self, type, reference):  # TODO alt. names, databases
        self.type = type
        self.reference = reference

    def as_version(self, version_number):
        return Version(self.type, self.reference, version_number)

    def __repr__(self):
        return '{}:{}'.format(self.type, self.reference)


class Version(UniqueIdentifier):
    def __init__(self, type, reference, version_number):  # TODO alt. names, databases
        super(Version, self).__init__(type, reference)
        self.version_number = version_number

    def __repr__(self):
        return '{}:{}.{}'.format(self.type, self.reference, self.version_number)


class DatabaseXRef(UniqueIdentifier):
    pass
    # /db_xref="GI:296146166"
    # /db_xref="GeneID:853656"
