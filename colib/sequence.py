from itertools import chain, repeat
import logging


class OverlapError(Exception):
    """
    This :class:`Exception` is raised when a mutation is applied to a position in a sequence that no longer
    exists because of a previous mutation.

    In *strict mode*, an :class:`OverlapError` is fired more frequently, such as when a deletion is applied to a range
    that has previously been modified by an insertion. In general, *strict mode* should be used when working
    with data from a single source of variant calls.
    """
    pass


class TranslationTable(object):
    """
    :param size: the length of the source sequence

    This class is inspired by the UCSC chain format for pairwise alignments documented here:

        http://genome.ucsc.edu/goldenPath/help/chain.html


    :class:`TranslationTable` encodes an alignment between two sequences, *source* and *target*, as a result
    of a sequence of mutations. Each mutation is either a `delete`, `insert`, or `substitute` operation where
    `substitute` is shorthand for a deletion followed by an insertion of equal size.

    .. attribute:: source_start

        The first position in the source sequence that aligns with the target sequence

    .. attribute:: source_end

        The last position in the source sequence that aligns with the target sequence.
    """

    def __init__(self, size):
        """
        test

        :param size: length of the source sequence
        """
        self.source_size = size
        self.target_size = size
        self.source_start, self.source_end = 0, size
        self.target_start, self.target_end = 0, size
        self.chain = [(size, 0, 0)]  # (ungapped_size, ds, dt)

    def invert(self):
        """
        Creates a copy of the :class:`TranslationTable` where *source* and *target* are inverted.

        :returns: a new :class:`TranslationTable` object.
        """
        tti = TranslationTable(self.target_size)
        tti.target_size = self.source_size
        tti.source_start, tti.source_end = self.target_start, self.target_end
        tti.target_start, tti.target_end = self.source_start, self.source_end
        tti.chain = [(ungapped_size, dt, ds) for ungapped_size, ds, dt in self.chain]
        return tti

    @classmethod
    def from_mutations(cls, sequence, mutations, strict=True):
        """
        :param sequence: the *source* sequence
        :param mutations: iterable of :class:`Mutation` objects
        :param strict: whether to use strict mode
        """
        tt = TranslationTable(len(sequence))

        for mutation in mutations:
            insert_size = len(mutation.new_sequence)
            translated_start = tt[mutation.start]

            if translated_start is None:
                raise OverlapError()

            if insert_size != mutation.size:
                tt.delete(translated_start, mutation.size, strict)
                if insert_size != 0:
                    tt.insert(translated_start, insert_size, strict)
            else:
                tt.substitute(translated_start, mutation.size, strict)

        return tt

    @classmethod
    def from_sequences(cls, reference, query, algorithm=None):
        raise NotImplementedError

    def _insert_gap(self, position, source_gap, target_gap, strict):
        if not 0 <= position <= self.source_size:
            raise IndexError()

        logging.debug('insert_gap: {} source gap: {}; target gap: {}'.format(position, source_gap, target_gap))

        if position <= self.source_start:
            first_ungapped_size, ds, dt = self.chain[0]
            if target_gap > first_ungapped_size:
                raise OverlapError('Cannot insert gap at start: overlaps')

            self.source_start += target_gap
            self.target_start += source_gap

            self.target_size += source_gap - target_gap
            # TODO decrease size of ungapped alignment

            self.chain[0] = first_ungapped_size - target_gap, ds, dt

            return

        elif position >= self.source_end:

            # only the exact position is legal for inserting at the end of the source sequence.
            if position == self.source_size:
                if target_gap:  # cannot delete what is not there:
                    raise IndexError()

                self.target_size += source_gap

                return

            # TODO OverlapError if non-strict

        else:
            offset = self.source_start

            for i, (ungapped_size, ds, dt) in enumerate(self.chain):
                print(i, position, ungapped_size, ds, dt)

                gap_start = offset + ungapped_size
                source_gap_end = gap_start + dt #- ds# NOTE might be ds

                #if position == offset:
                #    raise NotImplementedError()
                if offset < position < gap_start:
                    gap_offset = position - offset
                    ungapped_remainder = ungapped_size - gap_offset - target_gap

                    if ungapped_remainder < 0:
                        if strict:
                            raise OverlapError('Deletion at {} '
                                               'extends to following gap at {}'.format(position, position + target_gap))
                        else:
                            pass
                            # FIXME, eat up ungapped size up to the size of target_gap
                            #raise NotImplementedError()
                    elif ungapped_remainder == 0:  # merge chain tuples.

                        logging.debug(self.__dict__)
                        #raise NotImplementedError()

                    logging.warn(gap_offset)
                    self.chain.insert(i, (gap_offset, source_gap, target_gap))
                    self.chain[i + 1] = (ungapped_remainder, ds, dt)
                    break
                if position == gap_start:
                    # if ds and source_gap or dt and target_gap:
                    #     logging.debug(self.__dict__)
                    #     # allow for max one insertion and one deletion per coordinate:

                    # (source gap == insertion) an additional insertion is unusual and forbidden in strict mode
                    if strict and ds and source_gap:
                        raise OverlapError('Cannot insert gap at {}: '
                                           'Source sequence gap already present at this position.'.format(position))

                    # (target gap == deletion) deletion of an already deleted area is impossible and therefore raises
                    # an error
                    if dt and target_gap:
                        raise OverlapError('Cannot insert gap at {}: '
                                           'Target sequence gap already present at this position.'.format(position))

                    self.chain[i] = (ungapped_size, ds + source_gap, dt + target_gap)
                    break
                elif position < source_gap_end:
                    raise OverlapError('Cannot insert gap at {}: '
                                           'Gap already present from {} to {}'.format(position,
                                                                                      gap_start,
                                                                                      source_gap_end))

                offset = source_gap_end

            logging.debug('inserted gaps; source: %s, target: %s; target size: %s -> %s',
                          source_gap,
                          target_gap,
                          self.target_size,
                          self.target_size + source_gap - target_gap)

            # FIXME failure case first?
            self.target_size += source_gap - target_gap

            if source_gap:
                self.target_end += source_gap - target_gap
            if target_gap:
                # FIXME more complicated than this
                self.source_end += target_gap - source_gap

    def insert(self, position, size, strict=True):
        """
        When a sequence is inserted in the target sequence, this results in a gap in the alignment of the
        source sequence. The size of the target sequence also increases.

        The following edge cases occur:

        - If the insertion happens at position 0, no gap will be inserted, but the start of the source alignment
          is increased.
        - If the insertion happens at the end of the source sequence, likewise, the end of the source alignment
          increases. In both cases the size of the source sequence stays the same.
        - If the insertion happens inside a gap in the source sequence, or within the source sequence,
          the size of the source gap is increased unless `strict` is `True`, in which case an error is raised.
        - If the insertion happens inside a gap in the target sequence, an error is raised.


        """
        self._insert_gap(position, size, 0, strict)

    def delete(self, position, size, strict=True):
        self._insert_gap(position, 0, size, strict)

    def substitute(self, position, size, strict=True):
        self._insert_gap(position, size, size, strict)

    def le(self, position):
        """

        :func:`le()` attempts to return the coordinate in the target sequence that corresponds to the `position`
        parameter in the source sequence. If `position` falls into a gap in the target sequence, it will instead return
        the last coordinate in front of that gap.

        An `IndexError` is raised if the `position` does not exist in the source or if it maps to a coordinate before
        the start of the target sequence alignment.

        The result of the function is identical to the following implementation, but is more efficient::

            while True:
                target_position = self[position]
                if target_position is not None:
                    return target_position
                position -= 1

        :returns: the first position, equal or lower than `position` that exists in the query sequence.
        """
        while True:
            target_position = self[position]
            if target_position is not None:
                return target_position
            position -= 1

    def ge(self, position):
        """
        :returns: the first position, equal or greater than `position` that exists in the query sequence.
        :raises: IndexError
        """
        while True:
            query_position = self[position]
            if query_position is not None:
                return query_position
            position += 1

    @property
    def total_ungapped_size(self):
        """
        The total length of the alignment between source and target.
        """
        return sum(ungapped_size for ungapped_size, ds, dt in self.chain)

    def alignment(self):
        """
        Returns an iterator yielding tuples in the form ``(source, target)``.

        This function should only be used for debugging.

        :returns: an iterator over all coordinates in both the source and target sequence.
        """
        def source():
            offset = self.source_start
            for ungapped_size, ds, dt in self.chain:
                for i in range(ungapped_size):
                    yield offset + i, ' '
                for i in range(dt):
                    yield offset + ungapped_size + i, '*'
                for i in range(ds):
                    yield None
                offset += ungapped_size + dt

        def target():
            offset = self.target_start
            for ungapped_size, ds, dt in self.chain:
                for i in range(ungapped_size):
                    yield offset + i, ' '
                for i in range(ds):
                    yield offset + ungapped_size + i, '*'
                for i in range(dt):
                    yield None
                offset += ungapped_size + ds

        return zip(
            chain(repeat(None, self.source_start), source()),
            chain(repeat(None, self.target_start), target())
        )

    def alignment_str(self):
        s_str, t_str = '', ''
        logging.debug(self.__dict__)

        for s, t in self.alignment():
            s_str += ' {:3d}{}'.format(*s) if s is not None else '   - '
            t_str += ' {:3d}{}'.format(*t) if t is not None else '   - '

        return s_str + '\n' + t_str

    def __getitem__(self, position):
        """
        Maps from a source coordinate to a target coordinate.

        :returns: The new coordinate as an `int`, if it exists; `None` otherwise.
        :raises IndexError: If the coordinate does not exist in the original system
        """
        if not 0 <= position < self.source_size:
            raise IndexError()

        if position < self.source_start: # beginning of query was deleted.
            return None

        if position >= self.source_end:  # end of query was deleted.
            return None

        offset = self.source_start
        target_offset =  self.target_start
        for ungapped_size, ds, dt in self.chain:
            gap_start = offset + ungapped_size
            source_gap_end = gap_start + ds
            target_gap_end = gap_start + dt # NOTE might be ds

            # print 'chain:', (ungapped_size, ds, dt)
            # print (ungapped_size, ds, dt, target_offset, position, offset)
            # print position < offset + ungapped_size, position < offset + ungapped_size + dt
            if position < gap_start:  # target position is in this ungapped block
                # FIXME the <= is wrong, but necessary for some tests to work. The worst kind of wrong.
                logging.debug('{}: {} {} {} {}'.format(
                    position,
                    ungapped_size,
                    ds,
                    dt,
                    target_offset + (position - offset)
                ))
                return target_offset + (position - offset)
            elif position < target_gap_end:
                return None  # position falls into a gap in the query sequence (i.e. a deletion)
            # elif position < source_gap_end:
            else:
                target_offset += ungapped_size + ds
                offset += ungapped_size + dt
            # print

        logging.debug('offset, target_offset, position: ' + ', '.join(map(str, [offset, target_offset, position])))
        return target_offset + (position - offset)

    def __len__(self, other):
        """
        :returns: the length of the sequence in the source coordinate system.
        """
        return self.source_size
