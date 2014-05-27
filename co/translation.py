from itertools import chain, repeat
import logging
from Bio.SeqFeature import CompoundLocation, FeatureLocation


class OverlapError(Exception):
    """
    :class:`OverlapError` is raised when a mutation is applied to a position in a sequence that has been altered
    by a previous mutation.

    In *strict mode*, an :class:`OverlapError` is fired more frequently, such as when a deletion is applied to a range
    that has previously been modified by an insertion.
    """
    pass


class TranslationTable(object):
    """
    This class is inspired by the UCSC chain format for pairwise alignments documented here:

    http://genome.ucsc.edu/goldenPath/help/chain.html

    :class:`TranslationTable` encodes an alignment between two sequences, *source* and *target*.

    The alignment is encoded in a chain of tuples in the format ``(ungapped_size, ds, dt)``, where
    ``ungapped_size`` refers to regions that align, and the gaps ``dt`` and ``ds`` each refer to regions present
    only in the other sequence.

    .. attribute:: source_start

        The first position in the source sequence that aligns with the target sequence

    .. attribute:: source_end

        The last position in the source sequence that aligns with the target sequence.

    .. attribute: source_size

    .. attribute: target_start

    .. attribute: target_end

    .. attribute: target_size

    .. attribute: chain

    """
    def __init__(self,
                 source_size,
                 target_size,
                 source_start,
                 source_end,
                 target_start,
                 target_end,
                 chain):
        """
        test

        :param size: length of the source sequence
        """
        self.source_size = source_size
        self.target_size = target_size
        self.source_start, self.source_end = source_start, source_end
        self.target_start, self.target_end = target_start, target_end
        self.chain = chain

    def __invert__(self):
        return self.invert()

    def invert(self):
        """
        Creates a copy of the table where *source* and *target* are inverted.

        :returns: a new :class:`TranslationTable` object.
        """
        cls = self.__class__
        tti = cls.__new__(cls)
        tti.source_size = self.target_size
        tti.target_size = self.source_size
        tti.source_start, tti.source_end = self.target_start, self.target_end
        tti.target_start, tti.target_end = self.source_start, self.source_end
        tti.chain = [(ungapped_size, dt, ds) for ungapped_size, ds, dt in self.chain]
        return tti

    def le(self, position):
        """
        :func:`le()` attempts to return the coordinate in the target sequence that corresponds to the `position`
        parameter in the source sequence. If `position` falls into a gap in the target sequence, it will instead return
        the last coordinate in front of that gap.

        :raises IndexError: if `position` does not exist in the source or if it maps to a coordinate before
            the start of the target sequence alignment.
        :returns: the first position, equal or lower than `position` that exists in the query sequence.
        """
        while True:  # TODO a more efficient le() implementation
            target_position = self[position]
            if target_position is not None:
                return target_position
            position -= 1

    def ge(self, position):
        """
        .. seealso::

            The :func:`.le()` function.

        :raises IndexError: if `position` does not exist in the source or if it maps to a coordinate after
            the end of the target sequence alignment.
        :returns: the first position, equal or greater than `position` that exists in the query sequence.
        """
        while True:  # TODO a more efficient ge() implementation
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
            chain(repeat(None, self.target_start), target()))

    def alignment_str(self):
        """
        Returns a string representation of the alignment between `source` and `target` coordinates.

        .. warning::
            This function should only be used for debugging purposes.

        """
        s_str, t_str = '', ''
        for s, t in self.alignment():
            s_str += ' {:3d}{}'.format(*s) if s is not None else '   - '
            t_str += ' {:3d}{}'.format(*t) if t is not None else '   - '

        return s_str + '\n' + t_str

    def __iter__(self):
        for _ in range(self.source_start):
            yield None

        offset = self.source_start
        target_offset = self.target_start
        for ungapped_size, ds, dt in self.chain:
            for i in range(ungapped_size):
                yield target_offset + i
            for i in range(dt):
                yield None
            target_offset += ungapped_size + ds
            offset += ungapped_size + dt

        for _ in range(offset, self.source_size):
            yield None

    def __getitem__(self, position):
        """
        Maps from a source coordinate to a target coordinate.

        :returns: The new coordinate as an `int`, if it exists; `None` otherwise.
        :raises IndexError: If the coordinate does not exist in the original system
        """
        # TODO support slices (and optimize for them)

        if not 0 <= position < self.source_size:
            raise IndexError()

        if position < self.source_start:  # beginning of query was deleted.
            return None

        if position >= self.source_end:  # end of query was deleted.
            return None

        offset = self.source_start
        target_offset = self.target_start
        for ungapped_size, ds, dt in self.chain:
            gap_start = offset + ungapped_size
            source_gap_end = gap_start + ds
            target_gap_end = gap_start + dt  # NOTE might be ds

            if position < gap_start:  # target position is in this ungapped block
                # logging.debug('{}: {} {} {} {}'.format(
                #     position,
                #     ungapped_size,
                #     ds,
                #     dt,
                #     target_offset + (position - offset)
                # ))
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

    def __len__(self):
        """
        :returns: the length of the sequence in the source coordinate system.
        """
        return self.source_size


class MutableTranslationTable(TranslationTable):
    """
    :param size: the length of the source sequence

    A mutable version of :class:`TranslationTable` with `insert`, `delete` and `substitute` methods for updating
    the translation table with the corresponding mutations.

    """

    def __init__(self, size):
        """
        :param size: length of the source sequence
        """
        super(MutableTranslationTable, self).__init__(source_size=size,
                                                      target_size=size,
                                                      source_start=0,
                                                      source_end=size,
                                                      target_start=0,
                                                      target_end=size,
                                                      chain=[(size, 0, 0)])

    def freeze(self):
        """
        Return an immutable version of this translation table.

        :rtype: :class:`TranslationTable`
        """
        tt = TranslationTable.__new__(TranslationTable)
        tt.source_size = self.source_size
        tt.target_size = self.target_size
        tt.source_start, tt.source_end = self.source_start, self.source_end
        tt.target_start, tt.target_end = self.target_start, self.target_end
        tt.chain = list(self.chain)
        return tt

    @classmethod
    def from_mutations(cls, sequence, mutations, strict=True):
        """
        :param sequence: the *source* sequence
        :param mutations: iterable of :class:`Mutation` objects
        :param bool strict: use strict mode
        """
        tt = MutableTranslationTable(len(sequence))

        for mutation in mutations:
            insert_size = len(mutation.new_sequence)
            translated_start = tt[mutation.start]

            if translated_start is None:
                if strict:
                    raise OverlapError()

                translated_start = tt.le(translated_start)

            if insert_size != mutation.size:
                tt.delete(translated_start, mutation.size, strict)
                if insert_size != 0:
                    tt.insert(translated_start, insert_size, strict)
            else:
                tt.substitute(translated_start, mutation.size, strict)

        return tt

    @classmethod
    def from_sequences(cls, reference, query, algorithm=None):
        """
        :raises NotImplementedError:
        """
        raise NotImplementedError()

    def _insert_gap(self, position, source_gap, target_gap, strict):
        if not 0 <= position <= self.source_size:
            raise IndexError()
        if not position + target_gap - 1 <= self.source_size:
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
                gap_start = offset + ungapped_size
                source_gap_end = gap_start + dt  #- ds# NOTE might be ds

                logging.debug(
                    'source_gap_end={}; ({}) gap_start={}, dt={}, ds={}'.format(source_gap_end,
                                                                                position,
                                                                                gap_start,
                                                                                dt,
                                                                                ds))

                if offset < position < gap_start:
                    gap_offset = position - offset
                    ungapped_remainder = ungapped_size - gap_offset - target_gap

                    logging.debug('ungapped_remainder={}'.format(ungapped_remainder))

                    if ungapped_remainder < 0:
                        # TODO merge with position == gap_start case here.
                        raise OverlapError('Deletion at {} '
                                           'extends to following gap at {}'.format(position, position + target_gap))

                    elif ungapped_remainder == 0:  # reuse existing chain tuple
                        self.chain[i] = (gap_offset, ds + source_gap, dt + target_gap)
                        break

                    self.chain.insert(i, (gap_offset, source_gap, target_gap))
                    self.chain[i + 1] = (ungapped_remainder, ds, dt)
                    break
                elif position == gap_start:
                    # (source gap == insertion) an additional insertion is unusual and forbidden in strict mode
                    if strict and ds and source_gap:
                        raise OverlapError('Cannot insert gap at {}: '
                                           'Source sequence gap already present at this position.'.format(position))

                    # (target gap == deletion) deletion of an already deleted area is impossible and therefore raises
                    # an error
                    if dt and target_gap:
                        raise OverlapError('Cannot insert gap at {}: '
                                           'Target sequence gap already present at this position.'.format(position))

                    if target_gap:
                        if strict:
                            raise OverlapError('Strict gap insert at {} overlaps existing gap start.'.format(position))

                        # TODO change this code to an iterative version to support deletions that cross multiple gaps:
                        try:
                            next_ungapped_size, next_ds, next_dt = self.chain[i + 1]
                            next_ungapped_remainder = next_ungapped_size - target_gap

                            if next_ungapped_remainder == 0:
                                self.chain.pop(i + 1)
                                ds += next_ds
                                dt += next_dt
                            elif next_ungapped_remainder < 0:
                                raise OverlapError('Deletion at {} '
                                                   'extends through following gap at {}'.format(position,
                                                                                                position + target_gap))
                                # TODO make recursive here instead
                            else:
                                self.chain[i + 1] = (next_ungapped_size - target_gap, next_ds, next_dt)

                        except IndexError:
                            raise

                    self.chain[i] = (ungapped_size, ds + source_gap, dt + target_gap)

                    break
                elif position < source_gap_end:
                    logging.debug('position={}, {} {}; offset={}'.format(position, source_gap, target_gap, offset))
                    logging.debug(
                        'source_gap_end={}; gap_start={}, dt={}, ds={}'.format(source_gap_end, gap_start, dt, ds))
                    raise OverlapError('Cannot insert gap at {}: '
                                       'Gap already present from {} to {}'.format(position,
                                                                                  gap_start,
                                                                                  source_gap_end))
                elif position == source_gap_end:
                    next_ungapped_size, next_ds, next_dt = self.chain[i + 1]
                    ungapped_remainder = next_ungapped_size - target_gap

                    if ungapped_remainder < 0:
                        raise OverlapError('Deletion at {} '
                                           'extends to following gap at {}'.format(position, position + target_gap))
                    elif ungapped_remainder == 0:  # reuse existing chain tuple
                        self.chain[i] = (ungapped_size, ds + source_gap + next_ds, dt + target_gap + next_dt)
                        self.chain.pop(i + 1)
                    else:
                        self.chain[i] = (ungapped_size, ds + source_gap, dt + target_gap)
                        self.chain[i + 1] = (ungapped_remainder, next_ds, next_dt)
                    break

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
        Insert a gap in the *source* sequence

        :param int position: start of gap site
        :param int size: length of gap
        :param bool strict: use strict mode
        :raises OverlapError: in various edge cases involving overlapping mutations, particularly in `strict` mode.
        """
        self._insert_gap(position, size, 0, strict)

    def delete(self, position, size, strict=True):
        """
        Insert a gap in the *target* sequence

        :param int position: start of gap site
        :param int size: length of gap
        :param bool strict: use strict mode
        :raises OverlapError: in various edge cases involving overlapping mutations, particularly in `strict` mode.
        """
        self._insert_gap(position, 0, size, strict)

    def substitute(self, position, size, strict=True):
        """
        Insert two gaps of equal length in the *source* and *target* sequences

        :param int position: start of gap site
        :param int size: length of gap
        :param bool strict: use strict mode
        :raises OverlapError: in various edge cases involving overlapping mutations, particularly in `strict` mode.
        """
        self._insert_gap(position, size, size, strict)
