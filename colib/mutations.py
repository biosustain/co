from itertools import chain, repeat
import logging
import six
from colib.sequence import Sequence


class OverlapException(Exception):
    pass


class TranslationTable(object):
    """
    Inspired by the UCSC chain format for pairwise alignments.

    http://genome.ucsc.edu/goldenPath/help/chain.html
    """

    def __init__(self, size, changes=()):
        self._source_size = size
        self._target_size = size
        self._chain = [(size, 0, 0)]  # (ungapped_size, dr, dq)
        self._s_start, self._s_end = 0, size
        self._t_start, self._t_end = 0, size

    @classmethod
    def from_mutations(cls, sequence, mutations, strict=True):
        tt = TranslationTable(len(sequence))
        resolved_mutations = map(lambda m: m.resolve(sequence), mutations)

        for mutation in resolved_mutations:
            insert_size = len(mutation.new_sequence)
            translated_start = tt[mutation.start]

            if translated_start is None:
                raise OverlapException()

            if insert_size != mutation.size:
                tt.delete(translated_start, mutation.size)
                if insert_size != 0:
                    tt.insert(translated_start, insert_size)
            else:
                tt.substitute(translated_start, mutation.size)

        return tt

    @classmethod
    def from_sequences(cls, reference, query, algorithm=None):
        raise NotImplementedError

    def _insert_gap(self, position, source_gap, target_gap):
        if 0 > position >= self._source_size:
            raise IndexError()

        logging.debug('insert_gap: {} source gap: {}; target gap: {}'.format(position, source_gap, target_gap))

        if position <= self._s_start:
            first_ungapped_size, ds, dt = self._chain[0]
            if target_gap > first_ungapped_size:
                raise OverlapException('Cannot insert gap at start: overlaps')

            if source_gap:
                self._s_start += source_gap - target_gap
                self._t_end += source_gap - target_gap
            if target_gap:
                self._t_start += target_gap - source_gap
                self._s_end += target_gap - source_gap

            self._target_size += source_gap - target_gap
            # TODO decrease size of ungapped alignment

            self._chain[0] = first_ungapped_size - target_gap, ds, dt

            return

        elif position >= self._s_end:
            raise NotImplementedError()
        else:
            offset = self._s_start

            for i, (ungapped_size, ds, dt) in enumerate(self._chain):
                print(i, position, ungapped_size, ds, dt)

                gap_start = offset + ungapped_size
                source_gap_end = gap_start + dt # NOTE might be ds

                #if position == offset:
                #    raise NotImplementedError()
                if offset < position < gap_start:
                    gap_offset = position - offset
                    logging.warn(gap_offset)
                    self._chain.insert(i, (gap_offset, source_gap, target_gap))



                    self._chain[i + 1] = (ungapped_size - gap_offset , ds, dt)
                    break
                if position == gap_start:
                    if ds and source_gap or dt and target_gap:
                        # allow for max one insertion and one deletion per coordinate:
                        raise OverlapException('Cannot insert gap at {}: '
                                               'Gap of the same type already starting at this position.'.format(position))
                    self._chain[i] = (ungapped_size, ds + source_gap, dt + target_gap)
                    break
                elif position < source_gap_end:
                    raise OverlapException('Cannot insert gap at {}: '
                                           'Gap already present from {} to {}'.format(position,
                                                                                      gap_start,
                                                                                      source_gap_end))

                offset = source_gap_end

            # FIXME failure case first?
            self._target_size += source_gap - target_gap

            if source_gap:
                self._t_end += source_gap - target_gap
            if target_gap:
                self._s_end += target_gap - source_gap

    def insert(self, position, size):
        self._insert_gap(position, size, 0)

    def delete(self, position, size):
        self._insert_gap(position, 0, size)

    def substitute(self, position, size):
        self._insert_gap(position, size, size)

    @property
    def source_size(self):
        return self._source_size

    @property
    def target_size(self):
        return self._target_size

    def le(self, position): # FIXME horrible le implementation.
        """
        :returns: the first position, equal or lower than `position` that exists in the query sequence.
        """
        while True:
            query_position = self[position]
            if query_position is not None:
                return query_position
            position -= 1

    def ge(self, position): # FIXME horrible ge implementation.
        """
        :returns: the first position, equal or greater than `position` that exists in the query sequence.
        """
        while True:
            query_position = self[position]
            if query_position is not None:
                return query_position
            position += 1

    def alignment(self):
        """
        Returns an iterator yielding tuples in the form (source, target)

        :returns: an iterator over all coordinates in both the source and target sequence.
        """
        def source():
            offset = 0
            for ungapped_size, ds, dt in self._chain:
                for i in range(ungapped_size):
                    yield offset + i
                for i in range(dt):
                    yield offset + ungapped_size + i
                for i in range(ds):
                    yield None
                offset += ungapped_size + dt

        def target():
            offset = 0
            for ungapped_size, ds, dt in self._chain:
                for i in range(ungapped_size):
                    yield offset + i
                for i in range(ds):
                    yield offset + ungapped_size + i
                for i in range(dt):
                    yield None
                offset += ungapped_size + ds

        return zip(
            chain(repeat(None, self._s_start), source()),
            chain(repeat(None, self._t_start), target())
        )

    def alignment_str(self):
        s_str, t_str = '', ''

        for s, t in self.alignment():
            s_str += ' {:3d}'.format(s) if s is not None else ' ---'
            t_str += ' {:3d}'.format(t) if t is not None else ' ---'

        return s_str + '\n' + t_str

    def __getitem__(self, position):
        """
        Maps from a reference coordinate to a query coordinate.

        :returns: The new coordinate as an `int`, if it exists; none otherwise.
        :raises IndexError: If the coordinate does not exist in the original system
        """
        if position >= self._source_size or position < 0:
            raise IndexError()

        if position < self._t_start:  # beginning of query was deleted.
            return None

        if position >= self._t_end:  # end of query was deleted.
            return None

        offset = 0
        query_offset = self._s_start - self._t_start
        for ungapped_size, ds, dt in self._chain:
            gap_start = offset + ungapped_size
            source_gap_end = gap_start + ds
            query_gap_end = gap_start + dt # NOTE might be ds

            # print 'chain:', (ungapped_size, ds, dt)
            # print (ungapped_size, ds, dt, query_offset, position, offset)
            # print position < offset + ungapped_size, position < offset + ungapped_size + dt
            if position < gap_start:  # target position is in this ungapped block
                # FIXME the <= is wrong, but necessary for some tests to work. The worst kind of wrong.
                logging.debug('{}: {} {} {} {}'.format(
                    position,
                    ungapped_size,
                    ds,
                    dt,
                    query_offset + (position - offset)
                ))
                return query_offset + (position - offset)
            elif position < query_gap_end:
                return None  # position falls into a gap in the query sequence (i.e. a deletion)
            # elif position < source_gap_end:
            else:
                query_offset += ungapped_size + ds
                offset += ungapped_size + dt
            # print

        logging.debug('offset, query_offset, position: ' + ', '.join(map(str, [offset, query_offset, position])))
        return query_offset + (position - offset)

    def __len__(self, other):
        """
        Returns the length of the sequence in the old coordinate system.
        """
        return self._source_size


class Mutation(object):
    """

    DELINS does not have its own native Mutation class. These are best implemented in
    third-party libraries or using Mutation directly.

    A Mutation(start, end, new_sequence) is similar to the two derived mutations which however
    do not account for substitutions::

        DEL(start, end - start + 1), INS(start, new_sequence, replace=True)

    Mutation are stored as (position, size) pairs because (start, end) pairs do not allow for a proper zero-length
    mutation as it would be seen in an insertion. It is possible to simulate an insertion by keeping one character of
    the original sequence, but that would add ambiguity to the exact site of the mutated sequence.

    """
    def __init__(self, position, size, new_sequence=''):
        assert isinstance(new_sequence, (six.string_types, Sequence))
        self.position = int(position)
        self.size = int(size)
        self.new_sequence = new_sequence

    @property
    def start(self):
        return self.position

    @property
    def end(self): # FIXME *dangerous* needs review.
        if self.size in (0, 1):
            return self.position
        return self.position + self.size - 1

    @property
    def new_size(self):
        return len(self.new_sequence)

    def is_substitution(self):
        return self.size == len(self.new_sequence)

    def is_deletion(self):
        return self.size > len(self.new_sequence)

    def is_insertion(self):
        return not self.is_substitution() and self.new_size > 0

    def context(self):
        """
        The context of a mutation returns a description of possible feature annotations or parts that were used to
        create the mutation. The context is used to create human-readable feedback.
        """
        return dict(
            position=repr(self.position),
            size=repr(self.size),
            new_sequence=repr(self.new_sequence)  # FIXME sequence representation only if it is a Sequence object.
        )

    def __repr__(self):
        if self.is_insertion() and self.size == 0:
            return '<Mutation: at {} insert "{}">'.format(self.position, self.new_sequence)
        elif self.is_deletion() and self.new_size == 0:
            return '<Mutation: delete {}({})>'.format(self.position, self.size)
        else:
            return '<Mutation: change {}({}) to "{}">'.format(self.position, self.size, self.new_sequence)


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
    def __init__(self, pos, new_sequence, replace=False):
        """
        :param replace: if `True`, eliminates the original character at the position. Some variant call formats keep
            the first character of the original sequence in the replacement sequence.
        """
        super(INS, self).__init__(pos, int(replace), new_sequence)
