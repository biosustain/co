# coding: utf-8
import unittest

from colib import Component
from colib.mutations import DEL, INS


class FeatureTestCase(unittest.TestCase):
    def test_add_features(self):
        component = Component('GAGAGAGATATAGAGAGA')
        component.features.add(8, 4, name='tata', qualifiers={'a': 1})
        component.features.add(17, 1, name='end')
        component.features.add(0, 1, name='start')

        self.assertEqual(['start', 'tata', 'end'], [f.name for f in component.features])

    def test_inherit_features(self):
        component = Component('ABCDEFGHIerrorJKLMNOPQRSTUVXYZ')
        component.features.add(0, 3, name='abc')  # fine
        component.features.add(9, 5, name='error')  # fine
        component.features.add(6, 6, name='GHI..err')
        component.features.add(11, 6, name='ror..JKL')
        component.features.add(8, 7, name='I..error..J')  # fine
        component.features.add(29, 1, name='end')  # fine

        mutated = component.mutate([DEL(9, 5)])

        self.assertEqual({'JKL', 'Z', 'ABC', 'GHI', 'IJ'}, set(str(f.sequence) for f in mutated.features))

        # TODO insert mutation (need to test mutation.size = 0)
        # TODO snp mutation (need to test link breaking).

    def test_feature_non_strict_order(self):
        """
        In non-strict mode, potentially ambiguous mutations are allowed where the order in which
        the order in which the mutations are applied is significant.
        """
        component = Component('12345')
        component.features.add(2, 2, name='34')

        #   |████|      |████|
        # 12|3  4|5   12|3  4|5
        # 12|3xy4|5   12|3| - 5
        # 12|3xy|-5   12|3| xy5

        mutated_1 = component.mutate([DEL(3), INS(3, 'xy')], strict=False)
        self.assertEqual('3', str(list(mutated_1.features)[0].sequence))  # 3xy would also be acceptable.
        self.assertEqual('123xy5', str(mutated_1))

        mutated_2 = component.mutate([INS(3, 'xy'), DEL(3)], strict=False)
        self.assertEqual('3xy', str(list(mutated_2.features)[0].sequence))
        self.assertEqual('123xy5', str(mutated_2))

    def test_inherit_features_removed(self):
        pass

    def test_find_features(self):
        pass

    def test_feature_reverse_strand(self):
        pass