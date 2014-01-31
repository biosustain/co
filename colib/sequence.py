class Sequence(object):
    # FIXME delete this clss and instead use Bio.Seq.Seq

    @property
    def sequence(self):
        raise NotImplementedError

    @property
    def length(self):
        return len(self.sequence)

    # TODO slicers

    def __getitem__(self, item):
        return self.sequence[item]

    def __len__(self):
        return len(self.sequence)