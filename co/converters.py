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
            component.features \
                .add(feature.location,
                     type=feature.type,
                     strand=feature.strand,
                     id=feature.id,
                     qualifiers=feature.qualifiers,
                     ref=feature.ref,
                     ref_db=feature.ref_db)
        return component

    @classmethod
    def to_seq_record(cls, component):
        def convert_features(features):
            for feature in features:
                yield SeqFeature(feature.location,
                                 type=feature.type,
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
