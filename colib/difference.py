class Diff(object):

    def __init__(self, added, removed, changed=None):
        self.added = tuple(added)
        self.removed = tuple(removed)
        self.changed = tuple(changed) if changed else None