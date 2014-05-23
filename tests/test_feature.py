# coding: utf-8
import unittest
from Bio.SeqFeature import FeatureLocation

from colib import Component, Feature
from colib.mutation import DEL, INS


class FeatureTestCase(unittest.TestCase):
    def test_add_features(self):
        component = Component('GAGAGAGATATAGAGAGA')
        component.features.add(FeatureLocation(8, 12), qualifiers={'name': 'tata'})
        component.features.add(FeatureLocation(17, 18), qualifiers={'name': 'end'})
        component.features.add(FeatureLocation(0, 1), qualifiers={'name': 'start'})

        self.assertEqual(['start', 'tata', 'end'], [f.qualifiers.get('name') for f in component.features])

    def test_inherit_features(self):
        component = Component('ABCDEFGHIerrorJKLMNOPQRSTUVXYZ')
        component.features.add(FeatureLocation(0, 3), id='abc')  # fine
        component.features.add(FeatureLocation(9, 14), id='error')  # fine
        component.features.add(FeatureLocation(6, 12), id='GHI..err')
        component.features.add(FeatureLocation(11, 17), id='ror..JKL')
        component.features.add(FeatureLocation(8, 15), id='I..error..J')  # fine
        component.features.add(FeatureLocation(29, 30), id='end')  # fine

        mutated = component.mutate([DEL(9, 5)])

        self.assertEqual({'JKL', 'Z', 'ABC', 'GHI', 'IJ'}, set(str(f.seq) for f in mutated.features))
        # TODO insert mutation (need to test mutation.size = 0)
        # TODO snp mutation (need to test link breaking).

    def test_feature_non_strict_order(self):
        """
        In non-strict mode, potentially ambiguous mutations are allowed where the order in which
        the order in which the mutations are applied is significant.
        """
        component = Component('12345', feature_class=Feature)
        component.features.add(FeatureLocation(2, 4), type='34')

        #   |████|      |████|
        # 12|3  4|5   12|3  4|5
        # 12|3xy4|5   12|3  -|5
        # 12|3xy|-5   12|3 xy|5

        mutated_1 = component.mutate([DEL(3), INS(3, 'xy')], strict=False)
        self.assertEqual('3xy', str(list(mutated_1.features)[0].seq))  # 3 would also be acceptable.
        self.assertEqual('123xy5', str(mutated_1.seq))

        mutated_2 = component.mutate([INS(3, 'xy'), DEL(3)], strict=False)
        self.assertEqual('3xy', str(list(mutated_2.features)[0].seq))
        self.assertEqual('123xy5', str(mutated_2.seq))

    def test_inherit_features_removed(self):
        pass

    def test_find_features(self):
        pass

    def test_feature_reverse_strand(self):
        pass
