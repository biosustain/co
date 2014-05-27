class Diff(object):
    def __init__(self, added, removed, changed=None):
        self.added = tuple(added)
        self.removed = tuple(removed)
        self.changed = tuple(changed) if changed else None

    def __invert__(self):
        return Diff(added=self.removed, removed=self.added, changed=self.changed)

    def __len__(self):
        return len(self.added) + len(self.removed) + (not self.changed or len(self.changed))

    def __repr__(self):
        if self.changed:
            return 'Diff(added={}, removed={}, changed={})'.format(self.added, self.removed, self.changed)
        return 'Diff(added={}, removed={})'.format(self.added, self.removed)