import json
import re

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from co.component import Component
from co.feature import Feature
from co.identifiers import UniqueIdentifier


__all__ = (
    'Converter',
    'GenbankConverter',
    'FASTAConverter',
    'SBOLConverter',
    'JSONConverter'
)

GENBANK_META_ANNOTATIONS = ('accessions', 'comment', 'gi', 'organism', 'sequence_version', 'source', 'taxonomy')

# TODO look up SO:#-equivalent sequence ontology terms.
GENBANK_TYPE_SO_TERM_MAP = {
    '-10_signal': 'minus_10_signal',
    '-35_signal': 'minus_35_signal',
    '3_prime_UTR': 'three_prime_UTR',
    '5_prime_UTR': 'five_prime_UTR',
    'CAAT_signal': 'CAAT_signal',
    'CDS': 'CDS',
    'C_region': 'C_region',
    'D-loop': 'D_loop',
    'D_segment': None,
    'GC_signal': None,
    'J_segment': None,
    'LTR': 'LTR',
    'N_region': None,
    'RBS': 'RBS',
    'STS': 'STS',
    'S_region': None,
    'TATA_signal': None,
    'V_region': None,
    'V_segment': None,
    'assembly_gap': None,
    'attenuator': 'attenuator',
    'centromere': 'centromere',
    'enhancer': 'enhancer',
    'exon': 'exon',
    'gap': 'gap',
    'gene': 'gene',
    'iDNA': 'iDNA',
    'intron': 'intron',
    'mRNA': 'mRNA',
    'mat_peptide': 'mature_protein_region',
    'misc_RNA': 'RNA',
    'misc_binding': None,
    'misc_difference': None,
    'misc_feature': None,
    'misc_recomb': 'recombination_feature',
    'misc_signal': None,
    'misc_structure': None,
    'mobile_element': None,
    'modified_base': 'modified_DNA_base',
    'ncRNA': 'ncRNA',
    'old_sequence': None,
    'operon': 'operon',
    'oriT': 'oriT',
    'polyA_signal': 'polyA_signal_sequence',
    'polyA_site': 'polyA_site',
    'precursor_RNA': 'primary_transcript',
    'prim_transcript': 'primary_transcript',
    'primer_bind': 'primer_binding_site',
    'promoter': 'promoter',
    'protein_bind': 'protein_binding_site',
    'rRNA': 'rRNA',
    'rep_origin': 'origin_of_replication',
    'repeat_region': 'repeat_region',
    'sig_peptide': 'signal_peptide',
    'source': None,
    'stem_loop': 'stem_loop',
    'tRNA': 'tRNA',
    'telomere': 'telomere',
    'terminator': 'terminator',
    'tmRNA': 'tmRNA',
    'transit_peptide': 'transit_peptide',
    'unsure': None,
    'variation': None
}

SO_TERM_GENBANK_TYPE_MAP = {v: k for k, v in GENBANK_TYPE_SO_TERM_MAP.items()}

GENBANK_SINGLE_QUALIFIERS = (
    'gene',
    'locus_tag',
    'codon_start',
    'protein_id',
    'note',
    'pseudogene',
    'GO_component',
    'GO_function',
    'GO_process',
    'EC_number',
    'product',
    'experiment',
    'transl_table',
    'rpt_family',
)


class Converter(object):

    @classmethod
    def from_seq_record(cls, record):
        component = Component(
            seq=record.seq,
            parent=None,
            id=record.id,
            name=record.name,
            description=record.description,
            annotations=record.annotations)

        for feature in record.features:
            feature_qualifiers = feature.qualifiers

            # clean up synonyms for easier indexing:
            if 'gene_synonym' in feature.qualifiers:
                feature_qualifiers['gene_synonym'] = re.split(r';\s+?', feature.qualifiers['gene_synonym'][0])

            # flatten qualifiers where possible 'locus_id': ['b0001']  -> 'locus_id': 'b0001'
            for name in GENBANK_SINGLE_QUALIFIERS:
                if name in feature.qualifiers:
                    feature_qualifiers[name] = feature.qualifiers[name][0]

            component.features \
                .add(feature.location,
                     type=GENBANK_TYPE_SO_TERM_MAP.get(feature.type.replace("'", '_prime_')),
                     strand=feature.strand,
                     id=feature.id,
                     qualifiers=feature_qualifiers,
                     ref=feature.ref,
                     ref_db=feature.ref_db)
        return component

    @classmethod
    def to_seq_record(cls, component):
        def convert_features(features):
            for feature in features:
                yield SeqFeature(feature.location,
                                 type=SO_TERM_GENBANK_TYPE_MAP.get(feature.type, feature.type),
                                 strand=feature.strand,
                                 id=feature.id,
                                 qualifiers=feature.qualifiers,
                                 ref=feature.ref,
                                 ref_db=feature.ref_db)

        return SeqRecord(component.seq,
                         features=list(convert_features(component.features)),
                         id=component.id or "<unknown id>",
                         name=component.name or "<unknown name>",
                         description=component.description or "<unknown description>",
                         annotations=component.annotations)

    @classmethod
    def from_file(cls, file):
        raise NotImplementedError()

    @classmethod
    def to_file(cls, component, file):
        raise NotImplementedError()


class GenbankConverter(Converter):

    @classmethod
    def from_file(cls, file):
        record = SeqIO.read(file, 'genbank')
        component = Converter.from_seq_record(record)
        return component

    @classmethod
    def to_file(cls, component, file):
        return SeqIO.write(cls.to_seq_record(component), file, 'genbank')


class FASTAConverter(Converter):
    pass


class SBOLConverter(Converter):
    pass


class _ComponentEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Feature):
            return {
                'start': obj.start,
                'end': obj.end,
                'type': obj.type,
                'id': obj.id,
                'strand': obj.strand,
                'qualifiers': obj.qualifiers
            }
        elif isinstance(obj, Component):
            return {
                'parent': obj.parent.id if obj.parent else None,
                'sequence': str(obj.seq),
                'mutations': obj.diff(obj.parent) if obj.parent else None,
                'id': obj.id,
                'annotations': obj.annotations,
                'features': (tuple(obj.features.added), tuple(obj.features.removed))
            }
        elif isinstance(obj, UniqueIdentifier):
            return [obj.type, obj.reference]


class JSONConverter(Converter):
    @classmethod
    def to_file(cls, component, file, record_id=None):
        file.write(json.dumps(component, cls=_ComponentEncoder))
