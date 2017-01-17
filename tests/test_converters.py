from collections import Counter
import unittest

from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from co import Component
from co.converters import GenbankConverter


class GenbankConverterTestCase(unittest.TestCase):
    def test_phr332_import(self):
        filename = 'tests/fixtures/pbr322-sample-sequence.gb'

        component = GenbankConverter.from_file(filename)

        self.assertEqual(4196, len(component))
        self.assertEqual(19, len(component.features))
        self.assertEqual('TTCTCATGTT', str(component.seq[0:10]))

        self.assertEqual(
            Counter(
                {'-35_signal': 5,
                 '-10_signal': 5,
                 'CDS': 3,
                 'misc_RNA': 2,
                 'RBS': 2,
                 'rep_origin': 1,
                 'sig_peptide': 1}),
            Counter(feature.type for feature in component.features)
        )

        feature_0 = list(component.features)[0]
        feature_2 = list(component.features)[2]  # minus strand
        feature_15 = list(component.features)[15]  # minus strand

        self.assertEqual('TTGACA', str(feature_0.seq))
        self.assertEqual('TAAACT', str(feature_2.seq))
        self.assertEqual('ATGAGTATTCAACATTTCCGTGTCGCCC'
                         'TTATTCCCTTTTTTGCGGCATTTTGCCT'
                         'TCCTGTTTTTGCT', str(feature_15.seq))

    def test_export(self):
        component = Component(Seq('TATAGAGACACA', alphabet=DNAAlphabet()),
                              id='Magic_Brick',
                              annotations={'accession': 'MB1'})

        component.features.add(FeatureLocation(3, 7), type='-10_signal')
        component.features.add(FeatureLocation(7, 9), id='GA', type='CDS', qualifiers={'locus_id': 'b0001'})

        # GenbankConverter.to_file(component, 'fixtures/MB_1.gb', 'Magic_Brick')
        record = GenbankConverter.to_seq_record(component)
        self.assertEqual(open('tests/fixtures/MB_1.gb').read(), record.format('genbank'))

    maxDiff = None

    @unittest.SkipTest
    def test_export_combined(self):
        c1 = Component(Seq('AGAGAGAGAGA', alphabet=DNAAlphabet()), annotations={'organism': 'Predator'})
        c2 = Component(Seq('CGCGCGCCGCGCGCGCG', alphabet=DNAAlphabet()), annotations={'organism': 'Alien'})
        c3 = Component(Seq('AAAAAAAAAAATTTTAA', alphabet=DNAAlphabet()))
        c123 = Component.combine(c1, c2, c3)

        record = GenbankConverter.to_seq_record(c123)

        self.assertEqual(open('fixtures/c123.gb').read(), record.format('genbank'))