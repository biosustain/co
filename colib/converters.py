from Bio import SeqIO
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from colib.components import Component


GENBANK_META_ANNOTATIONS = ('accessions', 'comment', 'gi', 'organism', 'sequence_version', 'source', 'taxonomy')

# TODO look up SO:#-equivalent sequence ontology terms.
GENBANK_TYPE_SO_TERM_MAP = {
    '-10_signal': 'SO:0000175',
    '-35_signal': 'SO:0000176',
    '3_prime_UTR': 'SO:0000205',
    '5_prime_UTR': 'SO:0000204',
    'CAAT_signal': 'SO:0000172',
    'CDS': 'SO:0000316',
    'C_region': None,
    'D-loop': 'SO:0000297',
    'D_segment': None,
    'GC_signal': None,
    'J_segment': None,
    'LTR': 'SO:0000286',
    'N_region': None,
    'RBS': 'SO:0000552',
    'STS': 'SO:0000331',
    'S_region': None,
    'TATA_signal': None,
    'V_region': None,
    'V_segment': None,
    'assembly_gap': None,
    'attenuator': 'SO:0000140',
    'centromere': 'SO:0000577',
    'enhancer': 'SO:0000165',
    'exon': 'SO:0000147',
    'gap': 'SO:0000730',
    'gene': 'SO:0000704',
    'iDNA': 'SO:0000723',
    'intron': 'SO:0000188',
    'mRNA': 'SO:0000234',
    'mat_peptide': 'SO:0000419',
    'misc_RNA': 'SO:0000356',
    'misc_binding': None,
    'misc_difference': None,
    'misc_feature': None,
    'misc_recomb': 'SO:0000298',
    'misc_signal': None,
    'misc_structure': None,
    'mobile_element': None,
    'modified_base': 'SO:0000305',
    'ncRNA': 'SO:0000655',
    'old_sequence': None,
    'operon': 'SO:0000178',
    'oriT': 'SO:0000724',
    'polyA_signal': 'SO:0000551',
    'polyA_site': 'SO:0000553',
    'precursor_RNA': 'SO:0000185',
    'prim_transcript': 'SO:0000185',
    'primer_bind': 'SO:0005850',
    'promoter': 'SO:0000167',
    'protein_bind': 'SO:0000410',
    'rRNA': 'SO:0000252',
    'rep_origin': None,
    'repeat_region': 'SO:0000657',
    'sig_peptide': 'SO:0000418',
    'source': None,
    'stem_loop': 'SO:0000313',
    'tRNA': 'SO:0000253',
    'telomere': 'SO:0000624',
    'terminator': 'SO:0000141',
    'tmRNA': 'SO:0000584',
    'transit_peptide': 'SO:0000725',
    'unsure': None,
    'variation': None
}

SO_TERM_GENBANK_TYPE_MAP = {v: k for k, v in GENBANK_TYPE_SO_TERM_MAP.items()}

GENBANK_SINGLE_QUALIFIERS = (
    'locus_id',
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
                feature.qualifiers['gene_synonym'] = re.split(r';\s+?', feature.qualifiers['gene_synonym'][0])

            # flatten qualifiers where possible 'locus_id': ['b0001']  -> 'locus_id': 'b0001'
            for name in GENBANK_SINGLE_QUALIFIERS:
                if name in feature.qualifiers:
                    feature.qualifiers[name] = feature.qualifiers[name][0]

            component.features.add(feature.location.start,
                                   size=len(feature),
                                   type=GENBANK_TYPE_SO_TERM_MAP.get(feature.type),
                                   name=feature_name,
                                   qualifiers=feature_qualifiers)

    @classmethod
    def to_file(cls, component, file, record_id):
        record = cls.to_genbank_record(component, record_id)
        return SeqIO.write(record, file, 'genbank')

    @staticmethod
    def to_genbank_record(component, record_id):
        if not component.sequence:
            raise RuntimeError('Can only export records with an explicit sequence.')

        record = SeqRecord(component.sequence, id=record_id)
        record.annotations = component.qualifiers

        for feature in component.features:
            location = FeatureLocation(feature.start, feature.end)
            type = SO_TERM_GENBANK_TYPE_MAP.get(feature.type) or feature.type
            qualifiers = feature.qualifiers
            record.features.append(SeqFeature(location, type=type, qualifiers=qualifiers))
        return record


class SBOLConverter(object):

    @classmethod
    def from_file(self, file):
        raise NotImplementedError()

    @classmethod
    def to_file(self, component, file, record_id=None):
        raise NotImplementedError()
