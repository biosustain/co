class BasePosition(object):

    def value(self, component):
        """
        Resolves the position based on a parent component.

        :raises: IndexError if this position is not applicable to the component.
        """
        raise NotImplementedError()


class Position(BasePosition):

    def __init__(self, pos):
        self._pos = pos

    def value(self, component):
        return self._pos

    def __add__(self, other):
        """
        Syntactic sugar for creating offsets from position.

        :param other:
        :type other: `Distance` such as `_Feature.length`; or `int`.
        """
        return Offset(self, other)

    def __sub__(self, other):
        """
        Syntactic sugar for creating offsets or distances from positions.
        """


class Distance(object):
    def __init__(self, num):
        self._num = num


class Offset(BasePosition):
    """

    """
    DIRECTION_FORWARD, DIRECTION_REVERSE = 1, -1

    def __init__(self, position, offset, direction=DIRECTION_FORWARD):
        pass


class Range(object):

    def __init__(self, start, end):
        self._start = start
        self._end = end
        # assert start <= end