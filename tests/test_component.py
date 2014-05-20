import unittest
from Bio.Alphabet import Alphabet
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

import six

from colib import Component, Feature, Source
from colib.mutations import SNP, DEL, INS, Mutation, SUB


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

    def test_from_components_no_copy(self):
        a = Component('ABCDE', id='a')
        f = Component('FGH')
        i = Component('IJKLMNOPQ', id='i')

        combined = Component.combine(a, f, i)

        print(list(combined.features))

        self.assertEqual('ABCDEFGHIJKLMNOPQ', str(combined.seq))
        self.assertEqual([Feature(a, FeatureLocation(0, 5), ref='a'),
                          Feature(f, FeatureLocation(5, 8)),
                          Feature(i, FeatureLocation(8, 17), ref='i')], list(combined.features))

    def test_from_components_copy(self):
        co_1 = Component(Seq('AAATTTAAA'))
        co_1.features.add(FeatureLocation(3, 6), type='repeat', qualifiers={'name': 'ts'})
        co_2 = Component(Seq('G' * 10))
        co_2.features.add(FeatureLocation(0, 10), type='repeat', qualifiers={'name': 'gs'})

        combined = Component.combine(co_1, co_2, copy_features=True)

        self.assertEqual('AAATTTAAAGGGGGGGGGG', str(combined.seq))
        self.assertEqual([Feature(combined, FeatureLocation(3, 6), type='repeat', qualifiers={'name': 'ts'}),
                          Feature(combined, FeatureLocation(9, 19), type='repeat', qualifiers={'name': 'gs'})], list(combined.features))

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