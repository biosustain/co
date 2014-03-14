import six

from colib.sequence import Sequence


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
