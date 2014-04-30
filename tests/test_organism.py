import unittest
from Bio.Alphabet import Alphabet
import six
from colib import Component, Feature
from colib.converters import GenbankConverter, JSONConverter
from colib.mutations import SNP, Mutation
from colib.organisms import HaploidOrganism


class HaploidOrganismTestCase(unittest.TestCase):

    def test_list_features(self):

        c1 = Component('A' * 10)
        c1.features.add(1, 10, type='repeat')

        c2 = Component('T' * 10)
        c2.features.add(1, 10, type='repeat')

        self.assertEqual([Feature(c1, 1, 10, type='repeat')], list(c1.features))
        self.assertEqual([Feature(c2, 1, 10, type="repeat")], list(c2.features))

        s = HaploidOrganism(display_id='strain')
        s.add('a', c1)
        s.add('b', c2)
        s.add('c', c2)

        self.assertEqual(2, len(s.features))
        self.assertEqual([Feature(Component('AAAAAAAAAA', Alphabet()), 1, 10, type="repeat"),
                          Feature(Component('TTTTTTTTTT', Alphabet()), 1, 10, type="repeat")], list(s.features))

    def test_mutate_strain(self):
        strain = HaploidOrganism(display_id='Strain 1')

        feature_4 = strain.contigs[0].features[4]
        feature_5 = strain.contigs[0].features[5]

        replacement= Component('GAGAGA')

        mutations = [
            SNP(100, 'A'),
            SNP(105, 'G'),
            Mutation(40, 45, ''),
            Mutation(20, 21, 'GATGA'),
            Mutation(feature_4.start, end=feature_4.end + 20),
            Mutation(feature_5.start, end=feature_5.end, new_sequence=replacement)
        ]

        new_genome = strain.contigs[0].mutate(mutations)

        new_strain = HaploidOrganism('strain-1-1', parents=[strain])
        new_strain.contigs.add(new_genome) # replaces genome.

        self.storage.add_organism(new_strain) # saves strain.

        self.assertEqual(False, feature_4 in new_strain)
        self.assertEqual(True, feature_5 in new_strain)

        self.assertEqual(len(new_strain.diff(strain)), 2)

    @unittest.SkipTest
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
