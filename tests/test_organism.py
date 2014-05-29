import unittest

from Bio.Alphabet import Alphabet
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import six

from co import Component, Feature
from co.converters import GenbankConverter, JSONConverter
from co.difference import Diff
from co.mutation import SNP, Mutation
from co.organism import HaploidOrganism


class HaploidOrganismTestCase(unittest.TestCase):
    def test_list_features(self):
        c1 = Component('A' * 10)
        c1.features.add(FeatureLocation(1, 10), type='repeat')

        c2 = Component('T' * 10)
        c2.features.add(FeatureLocation(1, 10), type='repeat')

        self.assertEqual([Feature(c1, FeatureLocation(1, 10), type='repeat')], list(c1.features))
        self.assertEqual([Feature(c2, FeatureLocation(1, 10), type="repeat")], list(c2.features))

        strain = HaploidOrganism(id='strain')
        strain.set('a', c1)
        strain.set('b', c2)
        strain.set('c', c2)

        self.assertEqual(2, len(strain.features))
        self.assertEqual({Feature(Component('AAAAAAAAAA'), FeatureLocation(1, 10), type="repeat"),
                          Feature(Component('TTTTTTTTTT'), FeatureLocation(1, 10), type="repeat")}, set(strain.features))

    def test_mutate_strain(self):
        strain = HaploidOrganism(id='strain-1')

        genome = Component('A' * 25 + 'C' * 25 + 'G' * 25 + 'T' * 25)
        genome.features.add(FeatureLocation(0, 25), type='a')
        genome.features.add(FeatureLocation(25, 49), type='c')
        feature_3 = genome.features.add(FeatureLocation(50, 74), type='g')
        feature_4 = genome.features.add(FeatureLocation(75, 99), type='t')  # FIXME 75 should be possible, no?

        strain.set('genome', genome)

        mutations = [
            SNP(7, 'T'),
            SNP(30, 'G'),
            Mutation(31, 3, ''),
            Mutation(40, new_sequence='GATGA'),
            Mutation(feature_3.start, end=feature_3.end - 5),
            Mutation(feature_4.start, end=feature_4.end, new_sequence=Seq('GAGAGA'))
        ]

        new_genome = strain.components['genome'].mutate(mutations)

        new_strain = HaploidOrganism('strain-2', parent=strain)
        new_strain.set('genome', new_genome)  # sets genome.

        self.assertEqual(False, feature_3 in new_strain.features)
        self.assertEqual(False, feature_4 in new_strain.features)

        self.assertEqual(1, len(new_strain.diff(strain)))  # 'genome' changed
        self.assertEqual('AAAAAAATAAAAAAAAAAAAAAAAACCCCCGCCCCCCGATGACCCCCCCCCCGGGGGGAGAGA', str(new_genome.seq))

        genome_fdiff = new_strain.components['genome'].fdiff(strain.components['genome'])
        self.assertEqual(set(genome.features), set(genome_fdiff.removed))  # all features removed/changed
        self.assertEqual(Diff(changed=['genome']), new_strain.diff(strain))

        print(set(genome_fdiff.added))

        self.assertEqual({
                             Feature(new_genome, FeatureLocation(0, 25), type="a"),
                             Feature(new_genome, FeatureLocation(52, 56), type="g"),
                             Feature(new_genome, FeatureLocation(25, 50), type="c")}, set(genome_fdiff.added))

        self.assertEqual({('a', 'AAAAAAATAAAAAAAAAAAAAAAAA'),
                          ('c', 'CCCCCGCCCCCCGATGACCCCCCCC'),
                          ('g', 'GGGG')}, set((f.type, str(f.seq)) for f in genome_fdiff.added))

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
