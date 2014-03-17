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
             SNP(16, 'q'),                                   # ABCdEFG HIJKLMNOPqRSTU  VWXYZ
             DEL(1),                                         # A-CdEFG HIJKLMNOPqRSTU  VWXYZ
             INS(21, 'xx'),                                  # A-CdEFG HIJKLMNOPqRSTUxxVWXYZ
             Mutation(10, 9, 'oops'),                        # A-CdEFG HIJoops-----TUxxVWXYZ
             SUB(4, 'ef'),                                   # A-CdefG HIJoops-----TUxxVWXYZ
             Mutation(6, 1, 'Gg')],                          # A-CdefGgHIJoops-----TUxxVWXYZ
        strict=False)                                        # 0 1234567890123     457890123
                                                             # .    .    .         .   .


        self.assertEqual('ACdefGgHIJoopsTUxxVWXYZ', six.text_type(mutated.sequence))

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
        # FIXME self.assertEqual('01234599', six.text_type(Component('012345').mutate([INS(6, '99')]).sequence))

    @unittest.SkipTest
    def test_inherits(self):
        pass

    @unittest.SkipTest
    def test_diff(self):
        pass