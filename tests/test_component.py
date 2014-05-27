import unittest
from Bio.Alphabet import Alphabet
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

import six

from co import Component, Feature, Source
from co.mutation import SNP, DEL, INS, Mutation, SUB


class ComponentTestCase(unittest.TestCase):

    def test_mutate_1(self):
        component = Component('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        # .   .     .    .    .      .
        mutated = component.mutate(                             # 01234567 8901234567890  12345
                                     [SNP(3, 'd'),              # ABCdEFGH IJKLMNOPQRSTU  VWXYZ
                                      DEL(1),                   # A-CdEFGH IJKLMNOPQRSTU  VWXYZ
                                      INS(21, 'xx'),            # A-CdEFGH IJKLMNOPQRSTUxxVWXYZ
                                      Mutation(10, 9, 'oops'),  # A-CdEFGH IJoops-----TUxxVWXYZ
                                      SUB(4, 'ef'),             # A-CdefGH IJoops-----TUxxVWXYZ
                                     Mutation(7, 1, 'Hh')],     # A-CdefGHhIJoops-----TUxxVWXYZ
                                     strict=False)              # 0 1234567890123     457890123
                                                                # .    .    .         .   .

        self.assertEqual('ACdefGHhIJoopsTUxxVWXYZ', str(mutated.seq))

    def test_mutate_replace(self):
        self.assertEqual('01ttf2345', str(Component('012345').mutate([Mutation(1, 1, '1ttf')]).seq))
        self.assertEqual('0ott45', str(Component('012345').mutate([Mutation(1, 3, 'ott')]).seq))
        self.assertEqual('z12345', str(Component('012345').mutate([SNP(0, 'z')]).seq))

    def test_mutate_delete(self):
        self.assertEqual('01234', str(Component('012345').mutate([DEL(5)]).seq))
        self.assertEqual('01235', str(Component('012345').mutate([DEL(4)]).seq))
        self.assertEqual('0123', str(Component('012345').mutate([DEL(4, 2)]).seq))
        self.assertEqual('2345', str(Component('012345').mutate([DEL(0, 2)]).seq))

    def test_mutate_insert(self):
        self.assertEqual('99012345', str(Component('012345').mutate([INS(0, '99')]).seq))
        self.assertEqual('09912345', str(Component('012345').mutate([INS(1, '99')]).seq))
        self.assertEqual('9912345', str(Component('012345').mutate([INS(0, '99', replace=True)]).seq))
        # FIXME self.assertEqual('01234599', six.text_type(Component('012345').mutate([INS(6, '99')]).sequence))

    maxDiff = None

    def test_from_components_no_copy(self):
        a = Component('ABCDE', id='a')
        f = Component('FGH')
        i = Component('IJKLMNOPQ', id='i')

        combined = Component.combine(a, f, i)

        self.assertEqual('ABCDEFGHIJKLMNOPQ', str(combined.seq))

        qualifiers = {'mol_type': 'other DNA', 'organism': None}
        self.assertEqual([Feature(a, FeatureLocation(0, 5), type='source', qualifiers=qualifiers, ref='a'),
                          Feature(f, FeatureLocation(5, 8), type='source', qualifiers=qualifiers),
                          Feature(i, FeatureLocation(8, 17), type='source', qualifiers=qualifiers, ref='i')],
                         list(combined.features))

    def test_from_components_copy(self):
        co_1 = Component(Seq('AAATTTAAA'))
        co_1.features.add(FeatureLocation(3, 6), type='repeat', qualifiers={'name': 'ts'})
        co_2 = Component(Seq('G' * 10))
        co_2.features.add(FeatureLocation(0, 10), type='repeat', qualifiers={'name': 'gs'})

        combined = Component.combine(co_1, co_2, copy_features=True)

        self.assertEqual('AAATTTAAAGGGGGGGGGG', str(combined.seq))
        self.assertEqual([Feature(combined, FeatureLocation(3, 6), type='repeat', qualifiers={'name': 'ts'}),
                          Feature(combined, FeatureLocation(9, 19), type='repeat', qualifiers={'name': 'gs'})], list(combined.features))

    def test_lineage_simple(self):

        generation1 = Component('Ax')
        generation2 = generation1.mutate([INS(1, 'B')])
        generation3 = generation2.mutate([INS(2, 'C')])
        other = Component('a')

        self.assertEqual('ABx', str(generation2.seq))
        self.assertEqual(False, generation1.inherits_from(generation3))
        self.assertEqual(True, generation3.inherits_from(generation1))
        self.assertEqual(False, generation3.inherits_from(other))
        self.assertEqual([generation2, generation1], list(generation3.get_lineage()))


    def test_inherited_search(self):
        letters = Component('AABBDDEE', features=[
            SeqFeature(FeatureLocation(0, 1), type='vowel'),
            SeqFeature(FeatureLocation(2, 5), type='consonant'),
            SeqFeature(FeatureLocation(5, 6), type='vowel')])

        letters = letters.mutate([INS(4, 'CC')])

        self.assertEqual('AABBCCDDEE', str(letters.seq))
        self.assertEqual([Feature(letters, FeatureLocation(0, 1), type='vowel'),
                          Feature(letters, FeatureLocation(7, 8), type='vowel')], list(letters.features.find(type='vowel')))

        self.assertEqual([], list(letters.features.find(type='consonant', between_end=1)))

    def test_inherit_feature_delete(self):

        orig = Component('TTACCCATT', features=[SeqFeature(FeatureLocation(0, 1), type='tt')])
        f = orig.features.add(FeatureLocation(3, 6), type='ccc')

        self.assertEqual('CCC', str(f.seq))

        mutated = orig.mutate([DEL(3, 3)])

        self.assertEqual('TTAATT', str(mutated.seq))
        self.assertEqual([Feature(mutated, FeatureLocation(0, 1), type='tt')], list(mutated.features))

    def test_quickstart_feature_inherit(self):
        slogan = Component('CoPy is for DNA components', features=[
                         SeqFeature(FeatureLocation(0, 4), type='name'),
                         SeqFeature(FeatureLocation(12, 15), id='DNA')])

        self.assertEqual('components', str(slogan.features.add(FeatureLocation(16, 26)).seq))
        self.assertEqual(['CoPy', 'DNA', 'components'], [str(f.seq) for f in slogan.features])

        new_slogan = slogan.mutate([DEL(2, 2), DEL(12, 4)])

        self.assertEqual('Co is for components', str(new_slogan.seq))

        self.assertEqual([Feature(new_slogan, FeatureLocation(0, 2), type='name'),
                          Feature(new_slogan, FeatureLocation(10, 20))], list(new_slogan.features))
        self.assertEqual(['Co', 'components'], [str(f.seq) for f in new_slogan.features])


    @unittest.SkipTest
    def test_mutate_break_source(self):
        pass

    @unittest.SkipTest
    def test_inherits(self):
        pass

    @unittest.SkipTest
    def test_diff(self):
        pass

    @unittest.SkipTest
    def test_feature_diff(self):
        pass

    @unittest.SkipTest
    def test_eq(self):
        pass

    @unittest.SkipTest
    def test_get_lineage(self):
        pass