from collections import Counter
import unittest
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from colib import Component
from colib.converters import GenbankConverter


class GenbankConverterTestCase(unittest.TestCase):

    def test_phr332_import(self):
        filename = 'fixtures/pbr322-sample-sequence.gb'

        component = GenbankConverter.from_file(filename)

        self.assertEqual(4196, len(component))
        self.assertEqual(19, len(component.features))
        self.assertEqual('TTCTCATGTT', str(component[0:10]))

        self.assertEqual(
            Counter(
                {'minus_35_signal': 5,
                 'minus_10_signal': 5,
                 'CDS': 3,
                 'RNA': 2,
                 'RBS': 2,
                 'origin_of_replication': 1,
                 'signal_peptide': 1}),
            Counter(feature.type for feature in component.features)
        )

        feature_0 = list(component.features)[0]
        feature_2 = list(component.features)[2]    # minus strand
        feature_15 = list(component.features)[15]  # minus strand

        self.assertEqual('TTGACA', str(feature_0.sequence))
        self.assertEqual('TAAACT', str(feature_2.sequence))
        self.assertEqual('ATGAGTATTCAACATTTCCGTGTCGCCC'
                         'TTATTCCCTTTTTTGCGGCATTTTGCCT'
                         'TCCTGTTTTTGCT', str(feature_15.sequence))

    def test_export(self):
        component = Component('TATAGAGACACA', alphabet=DNAAlphabet(), meta={'accession': 'MB1'})
        component.features.add(3, 4, type='minus_10_signal')
        component.features.add(7, 2, name='GA', type='CDS', qualifiers={'locus_id': 'b0001'})

        # GenbankConverter.to_file(component, 'fixtures/MB_1.gb', 'Magic_Brick')
        record = GenbankConverter.to_genbank_record(component, 'Magic_Brick')
        self.assertEqual(open('fixtures/MB_1.gb').read(), record.format('genbank'))