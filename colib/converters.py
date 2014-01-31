from Bio import SeqIO


class GenbankConverter(object):

    def from_file(self, file):
        record = SeqIO.read(file, 'genbank')

    def to_file(self, component, file):
        pass


class SBOLConverter(object):

    def from_file(self, file):
        pass