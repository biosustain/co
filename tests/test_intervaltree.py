import unittest

from colib.interval import IntervalTree, Interval


class IntervalTreeTestCase(unittest.TestCase):
    def test_simple(self):
        it = IntervalTree()

        it.add(Interval(5, 10))
        it.add(Interval(-5, 0))
        it.add(Interval(1, 2))

        self.assertEqual(3, len(it))
        self.assertEqual([Interval(-5, 0), Interval(1, 2), Interval(5, 10)], list(it))

    def test_remove(self):
        it = IntervalTree()

        for i in range(5):
            it.add(Interval(i * 10, i * 10 + 5))

        it.remove(Interval(20, 25))
        it.remove(Interval(40, 45))

        self.assertEqual(3, len(it))
        self.assertEqual([Interval(0, 5), Interval(10, 15), Interval(30, 35)], list(it))

    def test_iter(self):
        it = IntervalTree()

        for i in range(5):
            it.add(Interval(i, i))

        for i, interval in enumerate(it):
            self.assertEqual(i, interval.start)

    def test_find_overlapping(self):
        it = IntervalTree()

        it.add(Interval(5, 10))
        it.add(Interval(15, 20))
        it.add(Interval(25, 30))
        it.add(Interval(35, 40))
        it.add(Interval(0, 40))
        it.add(Interval(10, 10))

        self.assertEqual([], list(it.find_overlapping(-10, -1)))
        self.assertEqual([Interval(0, 40)], list(it.find_overlapping(0, 1)))
        self.assertEqual([Interval(0, 40), Interval(5, 10)], list(it.find_overlapping(5, 6)))

        self.assertEqual([Interval(0, 40),
                          Interval(15, 20),
                          Interval(25, 30)], list(it.find_overlapping(20, 34)))

        self.assertEqual([Interval(0, 40),
                          Interval(5, 10),
                          Interval(10, 10)], list(it.find_overlapping(10, 10)))

    @unittest.SkipTest
    def test_size(self):
        pass