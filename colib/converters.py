import json
from Bio import SeqIO
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from colib.components import Component, _Feature
from colib.identifiers import UniqueIdentifier


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

class GenbankConverter(object):

    @classmethod
    def from_file(cls, file):
        record = SeqIO.read(file, 'genbank')

        record_meta = {name: record.annotations.get(name) for name in GENBANK_META_ANNOTATIONS}
        component = Component(record.seq, meta=record_meta)

        # TODO identifiers such as accession (with version) and gi number.

        for feature in record.features:
            if feature.type == 'source':
                continue

            feature_name = None
            feature_qualifiers = feature.qualifiers

            if 'gene' in feature.qualifiers:
                feature_name = feature.qualifiers['gene'][0]

            # clean up synonyms for easier indexing:
            if 'gene_synonym' in feature.qualifiers:
                feature_qualifiers['gene_synonym'] = re.split(r';\s+?', feature.qualifiers['gene_synonym'][0])

            # flatten qualifiers where possible 'locus_id': ['b0001']  -> 'locus_id': 'b0001'
            for name in GENBANK_SINGLE_QUALIFIERS:
                if name in feature.qualifiers:
                    feature_qualifiers[name] = feature.qualifiers[name][0]

            component.features\
                .add(feature.location.start,
                     size=len(feature),
                     strand=feature.location.strand,
                     type=GENBANK_TYPE_SO_TERM_MAP.get(feature.type.replace("'", '_prime_')),
                     name=feature_name,
                     qualifiers=feature_qualifiers)

        return component

    @classmethod
    def to_file(cls, component, file, record_id):
        record = cls.to_genbank_record(component, record_id)
        return SeqIO.write(record, file, 'genbank')

    @staticmethod
    def to_genbank_record(component, record_id):
        if not component.sequence:
            raise RuntimeError('Can only export records with an explicit sequence.')

        record = SeqRecord(component.sequence, id=record_id)
        record.annotations = component.meta

        for feature in component.features:
            location = FeatureLocation(feature.start, feature.end + 1, strand=feature.strand)
            type = SO_TERM_GENBANK_TYPE_MAP.get(feature.type) or feature.type
            qualifiers = feature.qualifiers
            record.features.append(SeqFeature(location, type=type.replace('_prime_', "'"), qualifiers=qualifiers))
        return record


class FASTAConverter(object):

    @classmethod
    def from_file(cls, file):
        raise NotImplementedError()

    @classmethod
    def to_file(cls, component, file, record_id=None):
        raise NotImplementedError()


class SBOLConverter(object):

    @classmethod
    def from_file(cls, file):
        raise NotImplementedError()

    @classmethod
    def to_file(cls, component, file, record_id=None):
        raise NotImplementedError()


class _ComponentEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, _Feature):
            return {
                'pos': obj.position,
                'size': obj.size,
                'type': obj.type,
                'name': obj.name,
                'strand': obj.strand,
                'qualifiers': obj.qualifiers,
                'link': [obj.link.id, obj.link_is_broken] if obj.link else None
            }
        elif isinstance(obj, Component):
            return {
                'parent': obj.parent.id if obj.parent else None,
                'sequence': str(obj.sequence),
                'mutations': obj.diff(obj.parent) if obj.parent else None,
                'id': obj.id,
                'meta': obj.meta,
                'features': obj.features.as_tuple()
            }
        elif isinstance(obj, UniqueIdentifier):
            return [obj.type, obj.identifier]

class JSONConverter(object):

    @classmethod
    def to_file(cls, component, file, record_id=None):
        file.write(json.dumps(component, cls=_ComponentEncoder))
