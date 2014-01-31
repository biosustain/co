import six
from colib.sequence import Sequence


class TranslationTable(object):
    """
    Inspired by the UCSC chain format for pairwise alignments.

    http://genome.ucsc.edu/goldenPath/help/chain.html
    """

    def __init__(self, size, changes=()):
        self._reference_size = size
        self._query_size = size
        self._chain = [(size, 0, 0)]  # (ungapped_size, dr, dq)
        self._r_start, self._r_end = 0, size
        self._q_start, self._q_end = 0, size

    @classmethod
    def from_mutations(cls, length, mutations, strict=True):
        pass

    @classmethod
    def from_sequences(cls, reference, query, algorithm=None):
        raise NotImplementedError

    def _insert_gap(self, position, reference_gap, query_gap):
        if 0 > position >= self._reference_size:
            raise IndexError()

        if position <= self._r_start:
            if reference_gap:
                self._r_start += reference_gap - query_gap
                self._q_end += query_gap - reference_gap
            if query_gap:
                self._q_start += query_gap - reference_gap
                self._r_end += query_gap - reference_gap

            self._query_size += reference_gap - query_gap
            # TODO increase size of ungapped alignment
            return

        offset = self._r_start
        for i, (ungapped_size, dr, dq) in enumerate(self._chain):
            if offset < position <= offset + ungapped_size:
                gap_offset = position - offset
                self._chain.insert(i, (gap_offset, reference_gap, query_gap))
                self._chain[i + 1] = (ungapped_size - gap_offset, dr, dq)

                # TODO change _r_end, _q_end
                return
            offset += ungapped_size + dr # NOTE might be dq

    def insert(self, position, size):
        self._insert_gap(position, size, 0)

    def delete(self, position, size):
        self._insert_gap(position, 0, size)

    def substitute(self, position, size):
        self._insert_gap(position, size, size)

    def __getitem__(self, position):
        """
        Maps from a reference coordinate to a query coordinate.

        :returns: The new coordinate as an `int`, if it exists; none otherwise.
        :raises IndexError: If the coordinate does not exist in the original system
        """
        if position >= self._reference_size or position < 0:
            raise IndexError()

        if position < self._q_start:  # beginning of query was deleted.
            return None
        #
        # if position >= self._q_end:  # beginning of query was deleted.
        #     return None

        offset = 0
        query_offset = self._r_start - self._q_start
        for ungapped_size, dr, dq in self._chain:
            # print 'chain:', (ungapped_size, dr, dq)
            # print (ungapped_size, dr, dq, query_offset, position, offset)
            # print position < offset + ungapped_size, position < offset + ungapped_size + dq
            if position < offset + ungapped_size:  # target position is in this ungapped block
                return query_offset + (position - offset)
            elif position < offset + ungapped_size + dq:
                return None  # position falls into a gap in the query sequence.

            query_offset += ungapped_size + dr
            offset += ungapped_size + dq
            # print





        return 0

    def __len__(self, other):
        """
        Returns the length of the sequence in the old coordinate system.
        """
        return self._reference_size


class Mutation(object):
    """

    DelIns does not have its own native Mutation class. These are best implemented in
    third-party libraries or using Mutation directly.

    """
    def __init__(self, start, end, new_sequence=''):
        assert start >= end
        assert isinstance(new_sequence, (six.string_types, Sequence))
        self._start = start
        self._end = end
        self._new_sequence = new_sequence

    @property
    def start(self):
        return int(self._start)

    @property
    def end(self):
        return int(self._end)

    @property
    def size(self):
        """
        .. note::

            This method is called size, not length, as it applies to the length of the deleted sequence for which there
            is no iterable.
        """
        return self._end - self._start + 1

    @property
    def new_sequence(self):
        return self._new_sequence

    def is_substitution(self):
        return self.size == len(self.new_sequence)

    def is_deletion(self):
        return self.size > len(self.new_sequence)

    def is_insertion(self):
        return len(self.new_sequence) > 0

    def context(self):
        """
        The context of a mutation returns a description of possible feature annotations or parts that were used to
        create the mutation. The context is used to create human-readable feedback.
        """
        return dict(
            start=repr(self.start),
            end=repr(self.end),
            new_sequence=repr(self.new_sequence)  # FIXME sequence representation only if it is a Sequence object.
        )


class SUB(Mutation):
    def __init__(self, pos, new_sequence):
        super(SUB, self).__init__(pos, pos + len(new_sequence) - 1, new_sequence)


class SNP(SUB):
    def __init__(self, pos, new_nucleotide):
        assert len(new_nucleotide) == 1
        super(SNP, self).__init__(pos, new_nucleotide)


class DEL(Mutation):
    def __init__(self, pos, size=1):
        super(DEL, self).__init__(pos, pos + size - 1)


class INS(Mutation):
    def __init__(self, pos, new_sequence):
        super(INS, self).__init__(pos, pos, new_sequence)
