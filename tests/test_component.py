import unittest
import six
from colib import Component
from colib.mutations import SNP, DEL, INS, Mutation, SUB


class ComponentTestCase(unittest.TestCase):

    def test_mutate_1(self):
        component = Component('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
                                                             # .   .     .    .    .      .
        mutated = component.mutate(                          # 0123456 78901234567890  12345
            [SNP(3, 'd'),                                    # ABCdEFG HIJKLMNOPQRSTU  VWXYZ
             DEL(1),                                         # A-CdEFG HIJKLMNOPQRSTU  VWXYZ
             INS(21, 'xx'),                                  # A-CdEFG HIJKLMNOPQRSTUxxVWXYZ
             Mutation(10, 9, 'oops'),                        # A-CdEFG HIJoops-----TUxxVWXYZ
             SUB(4, 'ef'),                                   # A-CdefG HIJoops-----TUxxVWXYZ
             Mutation(6, 1, 'Gg')],                          # A-CdefGgHIJoops-----TUxxVWXYZ
        strict=False)                                        # 0 1234567890123     457890123
                                                             # .    .    .         .   .

        self.assertEqual('ACdefGgHIJoopsTUxxVWXYZ', six.text_type(mutated.sequence))

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

    @unittest.SkipTest
    def test_inherits(self):
        pass

    @unittest.SkipTest
    def test_diff(self):
        pass