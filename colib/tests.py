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
import six
from colib.components import Component
from colib.organisms import UnicellularOrganism, Contig
from colib.storage import Storage
from colib.position import Range
from colib.mutations import SNP, Mutation, DEL, INS, TranslationTable, SUB


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

    @unittest.SkipTest
    def test_insert_end(self):
        pass # TODO

    @unittest.SkipTest
    def test_delete_end(self):
        pass # TODO

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

        mutated = component.mutate(                          # 0123456 78901234567890  1234
            [SNP(3, 'd'),                                    # ABCdEFG HIJKLMNOPQRSTU  VXYZ
             SNP(16, 'q'),                                   # ABCdEFG HIJKLMNOPqRSTU  VXYZ
             DEL(1),                                         # A-CdEFG HIJKLMNOPqRSTU  VXYZ
             INS(21, 'xx'),                                  # A-CdEFG HIJKLMNOPqRSTUxxVXYZ
             Mutation(10, 18, 'oops'),                       # A-CdEFG HIJoops-----TUxxVXYZ
             SUB(4, 'ef'),                                   # A-CdefG HIJoops-----TUxxVXYZ
             Mutation(6, 6, 'Gg')]                           # A-CdefGgHIJoops-----TUxxVXYZ
        )

        self.assertEqual('ACdefGgHIJoopsTUxxVXYZ', six.text_type(mutated.sequence))

    def test_mutate_replace(self):
        self.assertEqual('01ttf2345', six.text_type(Component('012345').mutate([Mutation(1, 1, '1ttf')]).sequence))
        self.assertEqual('0ott45', six.text_type(Component('012345').mutate([Mutation(1, 3, 'ott')]).sequence))
        self.assertEqual('z12345', six.text_type(Component('012345').mutate([SNP(0, 'z')]).sequence))

    def test_mutate_delete(self):
        self.assertEqual('01234', six.text_type(Component('012345').mutate([DEL(5)]).sequence))
        self.assertEqual('01235', six.text_type(Component('012345').mutate([DEL(4)]).sequence))
        self.assertEqual('0123', six.text_type(Component('012345').mutate([DEL(4, 2)]).sequence))
        self.assertEqual('2345', six.text_type(Component('012345').mutate([DEL(0, 2)]).sequence))

    def test_mutate_insert(self):
        self.assertEqual('99012345', six.text_type(Component('012345').mutate([INS(0, '99')]).sequence))
        self.assertEqual('09912345', six.text_type(Component('012345').mutate([INS(1, '99')]).sequence))
        self.assertEqual('9912345', six.text_type(Component('012345').mutate([INS(0, '99', replace=True)]).sequence))
        self.assertEqual('01234599', six.text_type(Component('012345').mutate([INS(6, '99')]).sequence))

    @unittest.SkipTest
    def test_inherits(self):
        pass

    @unittest.SkipTest
    def test_diff(self):
        pass


class FeatureTestCase(unittest.TestCase):

    def test_add_features(self):
        component = Component('GAGAGAGATATAGAGAGA')
        component.features.add(8, 4, name='tata', attributes={'a': 1})
        component.features.add(17, 1, name='end')
        component.features.add(0, 1, name='start')

        self.assertEqual(['start', 'tata', 'end'], [f.name for f in component.features])

    def test_inherit_features(self):
        component = Component('ABCDEFGHIerrorJKLMNOPQRSTUVXYZ')
        component.features.add(0, 3, name='abc')
        component.features.add(9, 5, name='error')
        component.features.add(6, 6, name='GHI..err')
        component.features.add(29, 1, name='end')

        mutated = component.mutate([DEL(9, 5)])

        print(mutated.sequence)

        for f in mutated.features:
            print(f)

        self.assertEqual(None, list(mutated.features))


    def test_inherit_features_removed(self):
        pass

    def test_find_features(self):
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


