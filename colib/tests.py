"""

Testing plan:

1. load genbank file
2. create strain
    2.1. apply breseq/variant call mutations
    2.2. apply feature removal, addition mutations
3. load SBOL part(s)
    3.1. apply part insertion mutations
        3.1.1. apply those relative to features
4. lineage
5. store/load from database
6. store/load components & mutations from JSON (for REST api)

"""
import unittest
from colib.components import Component
from colib.organisms import UnicellularOrganism, Contig
from colib.storage import Storage
from colib.position import Range
from colib.mutations import SNP, Mutation, DEL, INS, TranslationTable


class TranslationTableTestCase(unittest.TestCase):

    def setUp(self):
        self.sequence = 'ABCDEFGHIJ'
        self.tt = TranslationTable(len(self.sequence))

    def test_insert_start(self):
        # e.g. new sequence: 12345ABCDEFGHIJ
        self.tt.insert(0, 5)
        self.assertEqual(5, self.tt[0])
        self.assertEqual(6, self.tt[1])
        self.assertEqual(10, self.tt[5])
        self.assertEqual(14, self.tt[9])

        for i, expected_letter in enumerate(self.sequence):
            self.assertEqual(expected_letter, '12345ABCDEFGHIJ'[self.tt[i]])

    def test_delete_start(self):
        # e.g. new sequence: -----FGHIJ
        self.tt.delete(0, 5)
        self.assertEqual(None, self.tt[0])
        self.assertEqual(None, self.tt[4])
        self.assertEqual(0, self.tt[5])
        self.assertEqual(4, self.tt[9])

        for i, expected_letter in enumerate(self.sequence):
            if '-----FGHIJ'[i] == '-':
                self.assertIsNone(self.tt[i])
            else:
                self.assertEqual(expected_letter, 'FGHIJ'[self.tt[i]])

        with self.assertRaises(IndexError):
            _ = self.tt[10]

    def test_insert_end(self):
        pass

    def test_delete_end(self):
        pass

    def test_insert_mid(self):
        self.tt.insert(3, 2)
        # e.g.  ABCXXDEFGHIJ
        #       012  3456789
        #       012345678901
        self.assertEqual([0, 1, 2, 5, 6, 7, 8, 9, 10, 11], list(self.tt))

    def test_delete_mid(self):
        # e.g.     DEFG
        #       ABC----HIJ
        #       0123456789
        #       012    345
        self.tt.delete(3, 4)
        self.assertEqual([0, 1, 2, None, None, None, None, 3, 4, 5], list(self.tt))

    def test_multiple(self):
        pass

    def test_substitute(self):
        pass


class ComponentTestCase(unittest.TestCase):


    def test_mutate_1(self):

        component = Component('ABCDEFGHIJKLMNOPQRSTUVXYZ')

        mutated = component.mutate(                          # 0123456789012345678901234
            [SNP(3, 'd'),                                    # ABCdEFGHIJKLMNOPQRSTUVXYZ
             SNP(16, 'q'),                                   # ABCdEFGHIJKLMNOPqRSTUVXYZ
             DEL(1),                                         # A-CdEFGHIJKLMNOPqRSTUVXYZ
             INS(20, 'Uvx'),                                 # A-CdEFGHIJKLMNOPqRSTUvxYZ
             Mutation(Range(10, 18), 'oops')]                # A-CdEFGHIJoops-----TUvxYZ
        )

        self.assertEqual('ACdEFGHIJoopsTUvxYZ', mutated.sequence)

        pass

class StorageTestCase(unittest.TestCase):

    def setUp(self):
        self.storage = Storage()

    def tearDown(self):
        pass

    def test_add_component(self):
        component = Component()
        self.storage.components.add(component)


class UnicellularOrganismTestCase(unittest.TestCase):

    def setUp(self):
        self.storage = storage = Storage()
        storage.components.add(Component('genome-1'))
        storage.components.add(Component('part-1'))
        storage.components.add(Component('plasmid-1'))

        strain = UnicellularOrganism('strain-1')
        strain.contigs.add(Contig('genome', self.storage.components['genome-1']))
        strain.contigs.add(Contig('plasmid-a', self.storage.components['plasmid-1']))

        storage.add_organism(strain)

    def tearDown(self):
        pass

    def test_create_strain(self):
        strain = self.storage.organisms['strain-1']
        self.assertEqual(strain.as_dict(), {'display_id': 'strain-1', 'contigs': [{'name': 'genome', 'component': 'genome-1'}]})

    def test_mutate_strain(self):
        strain = self.storage.organisms['strain-1']
        feature_4 = strain.contigs[0].features[4]
        feature_5 = strain.contigs[0].features[5]

        mutations = [
            SNP(100, 'A'),
            SNP(105, 'G'),
            Mutation(Range(40, 45), None),
            Mutation(Range(20, 21), 'GATGA'),
            Mutation(Range(feature_4.start, feature_4.end + 20), None),
            Mutation(Range(feature_5.start, feature_5.start), self.storage.components['part-1'])
        ]

        new_genome = strain.contigs[0].mutate(mutations)

        new_strain = UnicellularOrganism('strain-1-1', parents=[strain])
        new_strain.contigs.add(new_genome) # replaces genome.

        self.storage.add_organism(new_strain) # saves strain.

        self.assertEqual(False, feature_4 in new_strain)
        self.assertEqual(True, feature_5 in new_strain)

        self.assertEqual(len(new_strain.diff(strain)), 2)


    def test_list_features(self):
        "Expect all features from genome-1 and plasmid-1."
        pass


