import unittest
from Bio.Alphabet import Alphabet

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

        self.assertEqual('ACdefGHhIJoopsTUxxVWXYZ', str(mutated))

    def test_mutate_replace(self):
        self.assertEqual('01ttf2345', str(Component('012345').mutate([Mutation(1, 1, '1ttf')])))
        self.assertEqual('0ott45', str(Component('012345').mutate([Mutation(1, 3, 'ott')])))
        self.assertEqual('z12345', str(Component('012345').mutate([SNP(0, 'z')])))

    def test_mutate_delete(self):
        self.assertEqual('01234', str(Component('012345').mutate([DEL(5)])))
        self.assertEqual('01235', str(Component('012345').mutate([DEL(4)])))
        self.assertEqual('0123', str(Component('012345').mutate([DEL(4, 2)])))
        self.assertEqual('2345', str(Component('012345').mutate([DEL(0, 2)])))

    def test_mutate_insert(self):
        self.assertEqual('99012345', str(Component('012345').mutate([INS(0, '99')])))
        self.assertEqual('09912345', str(Component('012345').mutate([INS(1, '99')])))
        self.assertEqual('9912345', str(Component('012345').mutate([INS(0, '99', replace=True)])))
        # FIXME self.assertEqual('01234599', six.text_type(Component('012345').mutate([INS(6, '99')]).sequence))

    def test_from_components_no_copy(self):
        a = Component('ABCDE')
        f = Component('FGH')
        i = Component('IJKLMNOPQ')

        combined = Component.from_components(a, f, i)

        self.assertEqual('ABCDEFGHIJKLMNOPQ', str(combined))
        self.assertEqual([Feature(combined, 0, 5, source=a),
                          Feature(combined, 5, 3, source=f),
                          Feature(combined, 8, 9, source=i)], list(combined.features))

    def test_from_components_copy(self):
        co_1 = Component('AAATTTAAA')
        co_1.features.add(3, 3, type='repeat', name='ts')
        co_2 = Component('G' * 10)
        co_2.features.add(0, 10, type='repeat', name='gs')

        combined = Component.from_components(co_1, co_2, copy_features=True)

        self.assertEqual('AAATTTAAAGGGGGGGGGG', str(combined))
        self.assertEqual([Feature(combined, 3, 3, type='repeat', name='ts'),
                          Feature(combined, 9, 10, type='repeat', name='gs')], list(combined.features))

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