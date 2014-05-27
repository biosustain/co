# coding: utf-8
from functools import reduce
import logging

from Bio import Alphabet
from Bio.Seq import Seq
import operator
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import six

from co.difference import Diff
from co.feature import Feature, ComponentFeatureSet, Source
from co.translation import OverlapError, MutableTranslationTable


class Component(object):
    """
    .. attribute:: features

        :class:`FeatureSet` containing the features present in this component.

    .. attribute:: id

        A unique identifier for this component; preferably either :class:`str` or :class:`UniqueIdentifier`.

    :param features: A list of additional features (features in addition to those inherited from ``parent``)
    :param removed_features: A list of removed features (features present in ``parent`` or one of its parents,
        but not present in this component)
    :param mutations: A list of mutations that have been applied to ``parent`` to arrive at ``seq``. The mutations will
        not be applied to ``seq`` again. Use :meth:`Component.mutate` to mutate a component.

    """
    def __init__(self,
                 seq,
                 parent=None,
                 features=None,
                 removed_features=None,
                 feature_class=None,
                 id=None,
                 name=None,
                 description=None,
                 annotations=None,
                 mutations=None):

        self.id = id
        self.name = name
        self.description = description

        # --- start of BioPython code ---
        # From SeqRecord.py, Biopython
        # Copyright 2000-2002 Andrew Dalke.
        # Copyright 2002-2004 Brad Chapman.
        # Copyright 2006-2010 by Peter Cock.
        # All rights reserved.
        # This code is part of the Biopython distribution and governed by its
        # license.  Please see the LICENSE file that should have been included
        # as part of this package.
        # annotations about the whole sequence
        if annotations is None:
            annotations = {}
        elif not isinstance(annotations, dict):
            raise TypeError("annotations argument should be a dict")
        self.annotations = annotations
        # --- end of BioPython code ---

        if isinstance(seq, six.string_types):
            seq = Seq(seq)

        self._seq = seq
        self._parent = parent
        self._mutations = mutations or ()
        self._mutations_tt = MutableTranslationTable.from_mutations(seq, self._mutations).freeze()
        self.features = ComponentFeatureSet(self, feature_class=feature_class)

        if features is not None:
            for feature in features:
                # FIXME need to translate SeqFeature to Feature!!!!!!!!!!
                self.features.add(self._convert_to_feature(feature))

        if removed_features is not None:
            for feature in removed_features:
                self.features.remove(feature)

    def __hash__(self):
        return hash(str(self.seq))  # Biopython uses id() to generate Seq.__hash__ and to test for equality

    def __len__(self):
        return len(self.seq)

    def __eq__(self, other):
        return isinstance(other, (Component, SeqRecord)) and str(self.seq) == str(other.seq)

    def __repr__(self):
        return '{}({}, id={})'.format(self.__class__.__name__, repr(self.seq), repr(self.id))

    def _convert_to_feature(self, feature):
        assert isinstance(feature, SeqFeature)
        if not isinstance(feature, Feature):
            return Feature(self,
                           type=feature.type,
                           location=feature.location,
                           location_operator=feature.location_operator,
                           strand=feature.strand,
                           id=feature.id,
                           qualifiers=feature.qualifiers,
                           ref=feature.ref,
                           ref_db=feature.ref_db)
        return feature

    @property
    def seq(self):
        """
        The nucleotide sequence of this component, type :class:`Bio.Seq`.
        """
        return self._seq

    def tt(self, ancestor=None):
        if ancestor in (None, self._parent):
            return self._mutations_tt
        raise KeyError("Translation table from {} to {} cannot be accessed.".format(repr(self), repr(ancestor)))

    @classmethod
    def combine(cls, *components, **kwargs):
        """
        Joins multiple components together, creating a source feature annotation for each.

        :param copy_features: whether to copy and translate all features from each of the components. If this feature
            is false; fresh ``source`` features will be created for each of the components.
        """
        copy_features = kwargs.pop('copy_features', False)

        if not components:
            return Component('')

        seq = reduce(operator.add, (c.seq for c in components))
        combined = Component(seq=seq)

        if copy_features:
            offset = 0
            for component in components:
                for feature in component.features:
                    combined.features.add(feature._shift(offset)) #move(offset + feature.position, component=combined))
                offset += len(component)
        else:
            offset = 0
            for component in components:
                organism = component.annotations.get('organism', component.annotations.get('source', None))
                combined.features.add(location=FeatureLocation(offset, offset + len(component.seq)),
                                      type='source',
                                      qualifiers={'organism': organism, 'mol_type': 'other DNA'},
                                      ref=component.id)

                offset += len(component)

        return combined


    @property
    def parent(self):
        return self._parent

    @staticmethod
    def _translate_feature(feature, component, tt=None):
        if tt is None:
            tt = component.tt(feature.component)
        #
        # logging.debug('TRANSLATE   {}'.format(feature))
        # logging.debug('TRANSLATE T  {}'.format(list(enumerate(tt))))

        offset = tt[feature.location.start] - feature.location.start

        # FIXME does not include ref and ref_db!
        # TODO create a "broken reference" object and re-attach here
        return feature._shift(offset, component)

    def mutate(self, mutations, strict=True):
        """
        Creates a copy of this :class:`Component` and applies all ``mutations`` in order.

        Mutations are applied in order based on the original coordinate system. If ``strict=True`` and multiple
        mutations affect the same area, a :class:`OverlapException` is raised. Otherwise, a :class:`OverlapException`
        may still be raised if a mutation falls into a position that has previously been deleted or the intended effect
        of a series of mutation is ambiguous.

        Although mutations are executed in order, the mutate function will attempt to resize features based on where
        their start or end based have been in the reference sequence if a base is lost and then restored over multiple
        mutations. For each deleted base, the algorithm tracks the next base to the left and to the right that was
        not deleted.

        :param strict: If ``True``, will fail if two mutations affect the same sequence.
        :param mutations: ``list`` of :class:`Mutation` objects
        :returns: A mutated :class:`Component`
        :raises OverlapException: When mutations overlap and interact in ambiguous ways that the algorithm can't
            handle.

        """
        component = Component(seq=self._seq, parent=self)  # TODO use __new__ instead
        features = component.features
        tt = MutableTranslationTable(size=len(self._seq))

        sequence = self.seq.tomutable()

        changed_features = set()

        # check that all mutations are in range:
        for mutation in mutations:
            if mutation.end > len(self.seq):
                raise IndexError('{} ends at {} but sequence length is {}'.format(mutation, mutation.end, len(self)))

        for mutation in mutations:
            # translate_to mutations, apply.
            # catch strict errors
            # flag features as edited/changed (copies are created as this happens), deleted.

            # FIXME handle INS at end of sequence.
            # FIXME requires proper handling of r/q_end to find out new size.

            logging.debug('')
            logging.debug('---> mutation: "%s"', mutation)

            if strict:
                translated_start = tt[mutation.position]
            else:
                try:
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

            affected_features = features.overlap(mutation.position, mutation.end)

            logging.debug('{} in {} to {} '.format(set(features), mutation.position, mutation.end))
            logging.debug('affected features: {}'.format(affected_features))

            # TODO also search previously affected features and exclude these.
            # TODO e.g. features = _FeatureSet(component, inherit=self.features)

            logging.debug('sequence: "%s"', sequence)

            if translated_start is None:
                raise OverlapError()

            logging.debug('new features: {}'.format(list(features)))
            logging.info('features in mutation: {}'.format(affected_features))

            for feature in affected_features | changed_features:
                # logging.debug(list(range(len(self.sequence))))
                # logging.debug(list(tt))
                logging.info('{} with sequence "{}" affected by {}.'.format(feature, feature.seq, mutation))

                if feature.end < mutation.start or feature.start > mutation.end:
                    continue

                try:
                    changed_features.remove(feature)
                except KeyError:
                    features.remove(feature)

                # TODO find untranslated equivalents and replace when all is over.
                if mutation.start > feature.start:
                    if mutation.end < feature.end - 1 or mutation.size == 0:  # mutation properly contained in feature.
                        logging.debug('FMMF from {} to {}({})'.format(
                            feature.start,  # tt.ge(feature.start),
                            feature.end,  # tt.le(feature.end),
                            feature.end - mutation.size + mutation.new_size))

                        changed_features.add(feature._move(feature.start,
                                                           feature.end - mutation.size + mutation.new_size))
                    else:
                        logging.debug('FMFM from {} to {}'.format(
                            tt[feature.start],
                            tt[mutation.start - 1]))  # tt[mutation.start] ?

                        changed_features.add(feature._move(feature.start, mutation.start))
                else:  # mutation.start >= feature.start

                    if mutation.end < feature.end - 1:
                        logging.debug('MFMF from {} to {}'.format(tt[mutation.end + 1], tt[feature.end]))

                        changed_features.add(feature._move(mutation.end + 1, feature.end))
                    else:
                        pass  # feature removed and not replaced.
                        logging.debug('MFFM feature removed')

            logging.debug('applying mutation: at %s("%s") translating from %s("%s") delete %s bases and insert "%s"',
                          translated_start,
                          sequence[translated_start],
                          mutation.start,
                          self.seq[mutation.start],
                          mutation.size,
                          mutation.new_sequence)

            # apply mutation to sequence:
            if mutation.new_size != mutation.size:  # insertion, deletion or delins
                del sequence[translated_start:translated_start + mutation.size]
                tt.delete(mutation.start, mutation.size, strict)

                if mutation.new_size != 0:  # insertion or DELINS
                    sequence = sequence[:translated_start] + mutation.new_sequence + sequence[translated_start:]
                    tt.insert(mutation.start, mutation.new_size, strict)
            else:  # substitution:
                if mutation.size == 1:  # SNP
                    sequence[translated_start] = mutation.new_sequence
                else:
                    sequence = sequence[:translated_start] \
                               + mutation.new_sequence \
                               + sequence[translated_start + mutation.new_size:]

                tt.substitute(mutation.start, mutation.size, strict)

            logging.debug('new sequence: "%s"', sequence)
            logging.debug('changed features: {}'.format(changed_features))

        logging.debug('mutated sequence: %s', str(sequence))

        component._seq = sequence.toseq()
        component._mutations = tuple(mutations)
        component._mutations_tt = tt
        component.features = features

        for feature in changed_features:
            logging.debug('changed: %s', feature)
            translated_feature = Component._translate_feature(feature, component, tt)
            features.add(translated_feature)

            logging.debug('translating %s: %s(%s) "%s" -> %s(%s) "%s"',
                          feature,
                          feature.start,
                          feature.end,
                          feature.seq,
                          translated_feature.start,
                          translated_feature.end,
                          translated_feature.seq)

        return component

    def get_lineage(self):
        """
        Iterate over the ancestors of this component.

        :return: iterator over :class:`Component` objects
        """
        component = self
        while component.parent:
            component = component.parent
            yield component

    def inherits_from(self, other):
        """
        Returns `True` if this object is a mutated version of `other`, `False` otherwise.

        :returns: ``bool``
        """
        return other in self.get_lineage()

    # noinspection PyProtectedMember
    def diff(self, other):
        """
        Returns a list of mutations that describe the difference between the sequence of this component and
        ``other``.
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
