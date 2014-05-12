import logging
import unittest

from colib import MutableTranslationTable
from colib.translation import OverlapError


class TranslationTableTestCase(unittest.TestCase):
    def setUp(self):
        self.sequence = 'ABCDEFGHIJ'
        self.tt = MutableTranslationTable(len(self.sequence))

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
        logging.debug(self.tt.__dict__)
        self.assertEqual([None, None, None, None, None, 0, 1, 2, 3, 4], list(self.tt))

        with self.assertRaises(IndexError):
            _ = self.tt[10]

    def test_insert_end(self):
        self.tt.insert(10, 5)

        logging.debug('\n' + self.tt.alignment_str())
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], list(self.tt))
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, None, None, None, None, None], list(self.tt.invert()))

        with self.assertRaises(IndexError):
            self.tt.insert(11, 1)

    def test_insert_before_end(self):
        self.tt.insert(9, 5)
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 14], list(self.tt))
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, None, None, None, None, None, 9], list(self.tt.invert()))
        self.assertEqual(15, self.tt.target_size)

    def test_delete_end(self):
        self.tt.delete(9, 1)

        logging.debug(self.tt.__dict__)
        self.assertEqual(9, self.tt.target_size)
        self.assertEqual(10, self.tt.source_size)
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, None], list(self.tt))
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8], list(self.tt.invert()))

        self.tt.delete(8, 1)
        self.assertEqual(8, self.tt.target_size)
        self.assertEqual(10, self.tt.source_size)
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, None, None], list(self.tt))
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7], list(self.tt.invert()))

    def test_delete_deleted(self):
        self.tt.delete(2, 5)
        self.assertEqual([0, 1, None, None, None, None, None, 2, 3, 4], list(self.tt))

        with self.assertRaises(OverlapError):
            self.tt.delete(3, 3)

        self.assertEqual([0, 1, None, None, None, None, None, 2, 3, 4], list(self.tt))

    def test_delete_next_deleted(self):
        logging.debug(self.tt.__dict__)
        self.tt.delete(2, 5)
        self.assertEqual([0, 1, None, None, None, None, None, 2, 3, 4], list(self.tt))

        logging.debug(self.tt.__dict__)
        logging.debug('\n' + self.tt.alignment_str())

        with self.assertRaises(OverlapError):
            self.tt.delete(2, 1)

        self.tt.delete(0, 1)

        #self.assertEqual([(1, 0, 6), (8, 0, 0)], self.tt._chain)
        #self.assertEqual([(1, 0, 1), (0, 0, 5), (8, 0, 0)], self.tt._chain)

        logging.debug(self.tt.__dict__)

        logging.debug('\n' + self.tt.alignment_str())
        self.assertEqual([None, 0, None, None, None, None, None, 1, 2, 3], list(self.tt))

    def test_delete_start_gap_delete(self):
        self.tt.delete(2, 5)
        self.assertEqual([0, 1, None, None, None, None, None, 2, 3, 4], list(self.tt))
        self.tt.delete(0, 1)
        logging.debug('\n' + self.tt.alignment_str())
        logging.debug(self.tt.__dict__)
        self.assertEqual([None, 0, None, None, None, None, None, 1, 2, 3], list(self.tt))

    def test_insert_mid(self):
        self.tt.insert(3, 2)

        # e.g.   ABC  DEFGHIJ
        # source 012  3456789
        # target 012345678901
        #        ABCXXDEFGHIJ

        self.assertEqual([0, 1, 2, 5, 6, 7, 8, 9, 10, 11], list(self.tt))

        self.assertEqual(12, self.tt.target_size)
        self.assertEqual(10, self.tt.total_ungapped_size)

        logging.debug(self.tt.__dict__)
        self.assertEqual([0, 1, 2, None, None, 3, 4, 5, 6, 7, 8, 9], list(self.tt.invert()))

    def test_start(self):
        pass

    def test_end(self):
        pass

    def test_ungapped_total(self):
        """
        `insert()`: generally the total ungapped size remains the same because the target sequence still aligns
        over the whole length of the source sequence.
        `delete()`: the total ungapped size shrinks because the target sequence no longer aligns to the source sequence.
        `substitute()`: like delete; the effect is similar to a deletion followed by an insertion (the
        substituted sequence no longer aligns)

        """
        pass

    def test_invert(self):
        pass

    def test_chain_merge_delete(self):
        pass

    def test_delete_mid(self):
        # e.g.     DEFG
        #        ABC----HIJ
        # source 0123456789
        # target 012    345

        self.tt.delete(3, 4)
        self.assertEqual([0, 1, 2, None, None, None, None, 3, 4, 5], list(self.tt))
        self.assertEqual([0, 1, 2, 7, 8, 9], list(self.tt.invert()))

        self.assertEqual(6, self.tt.target_size)
        self.assertEqual(6, self.tt.total_ungapped_size)

    def test_order_1(self):
        self.tt.delete(5, 1)  # ABCDE  -GHIJ
        self.assertEqual(9, self.tt.target_size)
        self.tt.insert(5, 2)  # ABCDExx-GHIJ
        self.assertEqual(11, self.tt.target_size)
        self.assertEqual([0, 1, 2, 3, 4, None, 7, 8, 9, 10], list(self.tt))

    def test_order_2(self):
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], list(self.tt))
        self.tt.insert(5, 2)  # ABCDExxFGHIJ
        self.assertEqual([0, 1, 2, 3, 4, 7, 8, 9, 10, 11], list(self.tt))
        self.assertEqual(self.tt.target_size, 12)
        self.tt.delete(5, 1)  # ABCDExx-GHIJ
        self.assertEqual(self.tt.target_size, 11)
        self.assertEqual([0, 1, 2, 3, 4, None, 7, 8, 9, 10], list(self.tt))

    def test_component_mutate_1(self):
        tt = MutableTranslationTable(26)
        self.assertEqual([(26, 0, 0)], tt.chain)
        self.assertEqual(26, tt.target_size)

        tt.substitute(3, 1)
        self.assertEqual([(3, 1, 1), (22, 0, 0)], tt.chain)

        tt.delete(1, 1)
        self.assertEqual([(1, 0, 1), (1, 1, 1), (22, 0, 0)], tt.chain)
        self.assertEqual(25, tt.target_size)

        tt.insert(21, 2)
        self.assertEqual([(1, 0, 1), (1, 1, 1), (17, 2, 0), (5, 0, 0)], tt.chain)

        tt.delete(10, 9)
        self.assertEqual([(1, 0, 1), (1, 1, 1), (6, 0, 9), (2, 2, 0), (5, 0, 0)], tt.chain)

        tt.insert(10, 4)
        self.assertEqual([(1, 0, 1), (1, 1, 1), (6, 4, 9), (2, 2, 0), (5, 0, 0)], tt.chain)

        tt.substitute(4, 2)
        self.assertEqual([(1, 0, 1), (1, 3, 3), (4, 4, 9), (2, 2, 0), (5, 0, 0)], tt.chain)

        tt.insert(6, 2)
        self.assertEqual([(1, 0, 1), (1, 5, 3), (4, 4, 9), (2, 2, 0), (5, 0, 0)], tt.chain)

        self.assertEqual([0, None, 1, None, None, None, 7, 8, 9, 10, None, None, None,
                          None, None, None, None, None, None, 15, 16, 19, 20, 21, 22, 23], list(tt))

    def test_zero_gap_merge(self):
        tt = MutableTranslationTable(10)

        tt.delete(5, 2)
        self.assertEqual([(5, 0, 2), (3, 0, 0)], tt.chain)
        self.assertEqual([0, 1, 2, 3, 4, None, None, 5, 6, 7], list(tt))

        tt.delete(3, 2)
        self.assertEqual([(3, 0, 4), (3, 0, 0)], tt.chain)
        self.assertEqual([0, 1, 2, None, None, None, None, 3, 4, 5], list(tt))
        self.assertEqual([0, 1, 2, 7, 8, 9], list(tt.invert()))

    def test_substitute(self):
        pass

    def test_le(self):
        tt = MutableTranslationTable(10)

        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [tt.le(i) for i in range(10)])

        tt.delete(5, 2)

        logging.debug(tt.__dict__)

        self.assertEqual([0, 1, 2, 3, 4, 4, 4, 5, 6, 7], [tt.le(i) for i in range(10)])

        tt.delete(9, 1)
        logging.debug(tt.__dict__)

        self.assertEqual([0, 1, 2, 3, 4, None, None, 5, 6, None], list(tt))
        self.assertEqual([0, 1, 2, 3, 4, 4, 4, 5, 6, 6], [tt.le(i) for i in range(10)])

        tt.delete(0, 1)

        self.assertEqual([None, 0, 1, 2, 3, None, None, 4, 5, None], list(tt))

        with self.assertRaises(IndexError):
            tt.le(0)

        self.assertEqual([0, 1, 2, 3, 3, 3, 4, 5, 5], [tt.le(i) for i in range(1, 10)])

    def test_ge(self):
        tt = MutableTranslationTable(10)

        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [tt.ge(i) for i in range(10)])

        tt.delete(1, 3)

        self.assertEqual([0, 1, 1, 1, 1, 2, 3, 4, 5, 6], [tt.ge(i) for i in range(10)])

        tt.delete(7, 3)

        for i in (7, 8, 9, 10):
            with self.assertRaises(IndexError):
                tt.ge(i)

        self.assertEqual([0, None, None, None, 1, 2, 3, None, None, None], list(tt))
        self.assertEqual([0, 1, 1, 1, 1, 2, 3], [tt.ge(i) for i in range(7)])