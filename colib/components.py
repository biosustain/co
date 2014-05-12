# coding: utf-8
from functools import reduce
import logging

from Bio import Alphabet
from Bio.Seq import Seq
import operator

from colib.difference import Diff
from colib.features import Feature, ComponentFeatureSet, Source
from colib.translation import OverlapError, MutableTranslationTable


class Component(Seq):
    """
    If the underlying library supports it, mutated `Component` objects can be stored as a set of the features
    that have been added or removed from the parent component. It is recommended that the library implements a
    caching strategy for quickly looking up the precise list of features contained in a component.

    `Component` objects, once created, SHOULDN'T be edited. However they MAY be edited until they are referred
    to by another component -- either by direct mutation or through a `_Feature`. A strategy for deleting components
    without destroying descendant objects may be necessary.

    .. attribute:: display_id

        A unique identifier for this component; preferably either :class:`str` or :class:`UniqueIdentifier`.

    .. seealso:: :class:`Bio.Seq`

    """

    def __init__(self,
                 sequence='',
                 alphabet=Alphabet.generic_alphabet,
                 parent=None,
                 meta=None,
                 display_id=None,
                 feature_class=Feature):

        if isinstance(sequence, Seq):
            sequence = str(sequence)

        super(Component, self).__init__(sequence, alphabet)

        # TODO support components without sequence.

        self._parent = parent
        self._sequence = sequence
        self._mutations = ()
        self._mutations_tt = MutableTranslationTable.from_mutations(self, self._mutations)
        self.display_id = display_id

        self.features = ComponentFeatureSet(self, feature_class=feature_class)
        self.meta = meta or {}

    def tt(self, ancestor=None):
        if ancestor in (None, self._parent):
            return self._mutations_tt
        raise KeyError("Translation table from {} to {} cannot be accessed.".format(repr(self), repr(ancestor)))

    @classmethod
    def from_components(cls, *components, **kwargs):
        """
        Joins multiple components together, creating a source feature annotation for each.

        .. warning::

            Optionally, features from the old components can be copied into the new component. This is strongly
            discouraged when working with large components and should only be done when assembling a component from
            temporary components that will later be discarded. A custom inspection can recursively enumerate features
            from linked components.

        :param copy_features: When `copy_features` is set, no features with links to the original components will be
            created.

        :param alphabet: Alphabet to use for the component. If not given, the alphabet for the first component
            will be used.
        """
        copy_features = kwargs.pop('copy_features', False)
        alphabet = kwargs.pop('alphabet', None)

        if not components:
            return Component('')

        sequence = reduce(operator.add, components)
        combined = Component(sequence=sequence, alphabet=alphabet or components[0].alphabet)

        if copy_features:
            offset = 0
            for component in components:
                for feature in component.features:
                    combined.features.add(feature.move(offset + feature.position, component=combined))
                offset += len(component)
        else:
            offset = 0
            for component in components:
                combined.features.add(offset, size=len(component), source=Source(component))
                offset += len(component)

        return combined


    @property
    def parent(self):
        return self._parent

    def clone(self):
        """
        Creates a copy of this :class:`Component`.

        The copy is linked to the same storage as the original, but is not saved automatically.
        """
        return Component(self._sequence, parent=self)

    def mutate(self, mutations, strict=True, clone=True):
        """
        Creates a copy of this :class:`Component` and applies all :param:`mutations`.

        Mutations are applied in order based on the original coordinate system. If `strict=True` and multiple
        mutations affect the same area, an `OverlapException` is raised. Otherwise, an `OverlapException` may
        still be raised if a mutation falls into a position that has previously been deleted.

        When a mutation affects a feature that has a linked component, that link will be flagged as `broken`.

        Memory
        ^^^^^^

        Although mutations are executed in order, the mutate function will attempt to resize features based on where
        their start or end based have been in the reference sequence if a based is lost and then restored in multiple
        mutations. For each deleted base, the algorithm tracks the next base to the left and to the right that was
        not deleted. In strict mode, since mutations are not allowed to find_overlapping, this only affects directly adjacent
        mutations.

        ::

              ████      ████          ██        ██     ██        ██          ██       ████
            123  45   123  45    123  45   123  45    123  45   123  45    12345    123  45
            123xy45   123  -5    123xy45   123  -5    123xy45   123  -5    123-5    123xx45
            123xy-5   123xy-5    123xy-5   123xy-5    123xy-5   123xy-5      █░       ████
              ███░      ███░          ░█        ░█     ██        ██

        :param strict: If `True`, will fail if two mutations affect the same sequence.
        :param clone: `True` by default; if `False`, the component will *NOT* be cloned and the changes will be applied
            to the original version instead.
        :returns: A mutated copy of the component (or a reference to the original component if `clone=False`)
        :raises OverlapException: if mutations find_overlapping
        """
        component = self.clone() if clone else self
        features = component.features
        sequence = self.tomutable()

        tt = MutableTranslationTable(len(self)) if clone else self._mutations_tt
        changed_features = set()

        logging.debug('original features: {}'.format(list(self.features)))
        logging.debug('new features: {}'.format(list(features)))

        # check that all mutations are in range:
        for mutation in mutations:
            if mutation.end >= len(self):
                raise IndexError('{} ends at {} but sequence length is {}'.format(mutation, mutation.end, len(self)))

        for mutation in mutations:
            # translate_to mutations, apply.
            # catch strict errors
            # flag features as edited/changed (copies are created as this happens), deleted.

            # FIXME handle addition to end of sequence.
            # FIXME requires proper handling of r/q_end to find out new size.

            logging.debug('')
            logging.debug('---> mutation: "%s"', mutation)

            if strict:
                translated_start = tt[mutation.position]
            else:
                try:
                    logging.debug(tt.__dict__)
                    translated_start = tt.ge(mutation.position)
                    logging.debug('translated start: nonstrict=%s; strict=%s', translated_start, tt[mutation.position])
                except IndexError:
                    try:
                        translated_start = tt.le(mutation.position)
                        logging.debug('translated start: nonstrict=%s; strict=%s', translated_start,
                                      tt[mutation.position])
                    except IndexError:
                        raise OverlapError("Cannot find non-strict mutation coordinate in mutated sequence.")

            # TODO strict mode implementation that also fires on other overlaps.

            # TODO features.
            affected_features = features.overlap(mutation.position, mutation.end)

            logging.debug('{} in {} to {} '.format(set(features), mutation.position, mutation.position + mutation.size))
            logging.debug('affected features: {}'.format(affected_features))

            # TODO also search previously affected features and exclude these.
            # TODO e.g. features = _FeatureSet(component, inherit=self.features)

            logging.debug('sequence: "%s"', sequence)
            logging.debug('alignment: \n%s', tt.alignment_str())

            if translated_start is None:
                raise OverlapError()

            logging.debug('new features: {}'.format(list(features)))
            logging.info('features in mutation: {}'.format(affected_features))

            for feature in affected_features | changed_features:
                # logging.debug(list(range(len(self.sequence))))
                # logging.debug(list(tt))
                logging.info('{} with sequence "{}" affected by {}.'.format(feature, feature.sequence, mutation))

                if feature.end < mutation.start or feature.start > mutation.end:
                    continue

                try:
                    changed_features.remove(feature)
                except KeyError:
                    features.remove(feature)

                # TODO find untranslated equivalents and replace when all is over.
                if mutation.start > feature.start:
                    if mutation.end < feature.end or mutation.size == 0:  # mutation properly contained in feature.
                        logging.debug('FMMF from {} to {}({})'.format(
                            feature.start,  # tt.ge(feature.start),
                            feature.end,  # tt.le(feature.end),
                            feature.end - mutation.size + mutation.new_size))

                        changed_features.add(feature.move(feature.start,
                                                          feature.size - mutation.size + mutation.new_size))
                    else:
                        logging.debug('FMFM from {} to {}'.format(
                            tt[feature.start],
                            tt[mutation.start - 1]))  # tt[mutation.start] ?

                        changed_features.add(feature.move(feature.start, mutation.start - feature.start))
                else:  # mutation.start >= feature.start

                    if mutation.end < feature.end:
                        logging.debug('MFMF from {} to {}'.format(tt[mutation.end + 1], tt[feature.end]))

                        changed_features.add(feature.move(mutation.end + 1, feature.end - mutation.end))
                    else:
                        pass  # feature removed and not replaced.
                        logging.debug('MFFM feature removed')

            logging.debug('applying mutation: at %s("%s") translating from %s("%s") delete %s bases and insert "%s"',
                          translated_start,
                          sequence[translated_start],
                          mutation.start,
                          self[mutation.start],
                          mutation.size,
                          mutation.new_sequence)

            # apply mutation to sequence:
            if mutation.new_size != mutation.size:  # insertion, deletion or delins
                del sequence[translated_start:translated_start + mutation.size]
                tt.delete(mutation.start, mutation.size)

                if mutation.new_size != 0:  # insertion or DELINS
                    sequence = sequence[:translated_start] + mutation.new_sequence + sequence[translated_start:]
                    tt.insert(mutation.start, mutation.new_size)
            else:  # substitution:
                if mutation.size == 1:  # SNP
                    sequence[translated_start] = mutation.new_sequence
                else:
                    sequence = sequence[:translated_start] \
                               + mutation.new_sequence \
                               + sequence[translated_start + mutation.new_size:]

                tt.substitute(mutation.start, mutation.size)

            logging.debug('new sequence: "%s"', sequence)
            logging.debug('changed features: {}'.format(changed_features))

        logging.debug('mutated sequence: %s', str(sequence))

        component._data = str(sequence.toseq())
        component._mutations = tuple(mutations)
        component._mutations_tt = tt
        component.features = features

        logging.debug('\n' + tt.alignment_str())

        for feature in changed_features:
            logging.debug('changed: %s', feature)
            translated_feature = feature.translate_to(component, tt)
            features.add(translated_feature)

            logging.debug('translating %s: %s(%s) "%s" -> %s(%s) "%s"',
                          feature,
                          feature.position,
                          feature.size,
                          feature.sequence,
                          translated_feature.position,
                          translated_feature.size,
                          translated_feature.sequence)

        return component

    def get_lineage(self):
        component = self
        while component.parent:
            component = component.parent
            yield component

    def inherits_from(self, other):
        """
        Returns `True` if this object is a mutated version of `other`, `False` otherwise.
        :returns: Boolean
        """
        return other in self.get_lineage()

    # noinspection PyProtectedMember
    def diff(self, other):
        """
        Returns the list of mutations that represent the difference between this component and another.

        .. note::

            This may be aided by the storage library. The storage library (if attached) should be called first to see if
            there is some sort of optimized lookup. The library itself will either implement an optimized version or
            fall back to the diff solution used here.
        """
        assert isinstance(other, Component)
        if other == self.parent:
            return self._mutations
        raise NotImplementedError()

    def fdiff(self, other):
        if other == self.parent:
            return Diff(added=self.features.added, removed=self.features.removed)
        elif other.parent == self:
            return ~other.fdiff(self)
        raise NotImplementedError()
