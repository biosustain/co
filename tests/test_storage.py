import unittest
import six
from colib import Component
from colib.converters import GenbankConverter, JSONConverter
from colib.organisms import HaploidOrganism
from colib.storage import Storage


class StorageTestCase(unittest.TestCase):

    def setUp(self):
        self.storage = Storage()

    def tearDown(self):
        pass

    def test_add_component(self):
        component = Component()
        self.storage.components.add(component)

    def test_yeast(self):
        yeast = HaploidOrganism('Saccharomyces cerevisiae S288c')

        components = (
            ('I', 'NC_001133.9.gb', 'chromosome'),
            ('II', 'NC_001134.8.gb', 'chromosome'),
            ('III', 'NC_001135.5.gb', 'chromosome'),
            ('IV', 'NC_001136.10.gb', 'chromosome'),
            ('V', 'NC_001137.3.gb', 'chromosome'),
            ('VI', 'NC_001138.5.gb', 'chromosome'),
            ('VII', 'NC_001139.9.gb', 'chromosome'),
            ('VIII', 'NC_001140.6.gb', 'chromosome'),
            ('IX', 'NC_001141.2.gb', 'chromosome'),
            ('X', 'NC_001142.9.gb', 'chromosome'),
            ('XI', 'NC_001143.9.gb', 'chromosome'),
            ('XII', 'NC_001144.5.gb', 'chromosome'),
            ('XIII', 'NC_001145.3.gb', 'chromosome'),
            ('XIV', 'NC_001146.8.gb', 'chromosome'),
            ('XV', 'NC_001147.6.gb', 'chromosome'),
            ('XVI', 'NC_001148.4.gb', 'chromosome'),
            ('MT', 'NC_001224.1.gb', 'mitochondrial_chromosome')
        )

        for name, path, type in components:
            yeast.components.add((name, GenbankConverter.from_file('fixtures/S288C/' + path), type))

        output = six.StringIO()
        JSONConverter.to_file(GenbankConverter.from_file('fixtures/S288C/NC_001133.9.gb'), output)

        self.assertEqual('', output.getvalue())

        self.storage.organisms.add(yeast)