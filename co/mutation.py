import six
from Bio.Seq import Seq

from co import Component


class Mutation(object):
    """
    A ``Mutation(start, size, new_sequence)`` is similar to the two derived mutations::

        DEL(start, size), INS(start, new_sequence)

    :param int position: start index
    :param int size: length of deletion
    :param new_sequence: insertion sequence
    :type new_sequence: :class:`str`, :class:`Component` or :class:`Bio.Seq`

    .. note::
        Mutation are stored as ``(position, size)`` pairs because ``(start, end)`` pairs do not allow for
        unambiguous zero-length mutations (i.e. insertions). It is possible to simulate an insertion by keeping one
        character of the original sequence, but that would introduce ambiguity to the exact site of the mutated sequence.

    .. attribute:: new_sequence

        Replacement sequence inserted at :attr:`position`.

    .. attribute:: size

        Length of the stretch of original sequence that is deleted at :attr:`position`.

    .. attribute:: position

        Start index of the mutation, zero-based.

    """

    def __init__(self, position, size=None, new_sequence='', end=None):
        assert isinstance(new_sequence, six.string_types + (Seq, Component))
        assert size is None or size >= 0
        assert end is None or end >= position

        self.position = int(position)

        if end is not None:
            self.size = int(end) - int(position) + 1
        elif size is None:
            self.size = 0
        else:
            self.size = int(size)

        if isinstance(new_sequence, Seq):
            new_sequence = str(new_sequence)

        self.new_sequence = new_sequence

    @property
    def start(self):
        """
        Identical to :attr:`Mutation.position`
        """
        return self.position

    @property
    def end(self):  # FIXME *dangerous* needs review.
        """
        Computed end coordinate of the deletion. Use with caution.
        """
        if self.size in (0, 1):
            return self.position
        return self.position + self.size - 1

    @property
    def new_size(self):
        """
        The length of :attr:`new_sequence`
        """
        return len(self.new_sequence)

    def is_substitution(self):
        """
        ``True`` if the size of the deletion is equal to the size of the insertion.
        """
        return self.size == len(self.new_sequence)

    def is_deletion(self):
        """
        ``True`` if the size of the deletion is larger than the size of the insertion
        """
        return self.size > len(self.new_sequence)

    def is_insertion(self):
        """
        ``True`` if the size of the deletion is zero and ``new_sequence`` is not empty.
        """
        return not self.is_substitution() and self.new_size > 0

    def __repr__(self):
        if self.is_insertion() and self.size == 0:
            return 'INS({}, "{}")'.format(self.position, self.new_sequence)
        elif self.is_deletion() and self.new_size == 0:
            return 'DEL({}, {})'.format(self.position, self.size)
        else:
            return 'Mutation({}, {}, "{}")'.format(self.position, self.size, self.new_sequence)


class SUB(Mutation):
    def __init__(self, pos, new_sequence):
        super(SUB, self).__init__(pos, len(new_sequence), new_sequence)


class SNP(SUB):
    def __init__(self, pos, new_nucleotide):
        assert len(new_nucleotide) == 1
        super(SNP, self).__init__(pos, new_nucleotide)


class DEL(Mutation):
    def __init__(self, pos, size=1):
        super(DEL, self).__init__(pos, size)


class INS(Mutation):
    """

    :param int pos: zero-based insertion index
    :param new_sequence: insertion sequence
    :type new_sequence: str or Bio.Seq
    :param bool replace: if ``True``, eliminates the original character at the position. Some variant call formats keep
        the first character of the original sequence in the replacement sequence.
    """

    def __init__(self, pos, new_sequence, replace=False):
        super(INS, self).__init__(pos, int(replace), new_sequence)
