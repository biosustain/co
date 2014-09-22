from abc import ABCMeta, abstractproperty

import six


__all__ = ('IntervalTree', 'BaseInterval', 'Interval')


def node_search(node, data):
    if node is None:
        return

    if node.left is not None:
        for d in node_search(node.left, data):
            yield d

    if node.data.overlaps(data):
        yield node.data

    # if node data > search data then the right tree will also be > search data
    if data.end < node.data.start:
        return

    if node.right is not None:
        for d in node_search(node.right, data):
            yield d


def node_traverse(node):
    if node is None:
        return

    for i in node_traverse(node.left):
        yield i

    yield node.data

    for i in node_traverse(node.right):
        yield i


class IntervalTree(object):
    """

    IntervalMixin Tree based on a binary search tree.

    """

    def __init__(self):
        self.root = None
        self.size = 0

    def add(self, interval):
        if self.root is None:
            self.root = IntervalNode(interval)
        else:
            current = self.root
            while True:
                if interval < current.data:
                    if current.left:
                        current = current.left
                    else:
                        current.left = IntervalNode(interval)
                        break
                else:
                    if current.right:
                        current = current.right
                    else:
                        current.right = IntervalNode(interval)
                        break

        self.size += 1  # FIXME since original node might be replaced, size will not necessarily change.

    def copy(self):
        it = IntervalTree()
        for node in self:
            it.add(it.data)
        return it

    def remove(self, interval):
        try:
            if self.root:
                self.root.delete(interval)
                self.size -= 1
        except LookupError:
            pass

    def find_overlapping(self, start, end):
        assert end >= start
        for i in node_search(self.root, Interval(start, end)):
            yield i

    def __iter__(self):
        stack = []
        node = self.root
        while stack or node:
            if node:
                stack.append(node)
                node = node.left
            else:  # we are returning so we pop the node and we yield it
                node = stack.pop()
                yield node.data
                node = node.right

    def __len__(self):
        return self.size


class BaseInterval(six.with_metaclass(ABCMeta, object)):
    """

    An interval ``Interval(start, end)`` describes a range `[start, end]` where both ``start`` and ``end`` are included.

    .. attribute:: start

    .. attribute:: end
    """

    @abstractproperty
    def start(self):
        raise NotImplementedError()

    @abstractproperty
    def end(self):
        raise NotImplementedError()

    def contains(self, other):
        return self.start <= other.start and self.end >= other.end

    def overlaps(self, other):
        return self.start <= other.end and self.end >= other.start

    def __contains__(self, position):
        return self.start <= position <= self.end

    def __lt__(self, other):
        if self.start < other.start:
            return True
        if self.start > other.start:
            return False
        if self.end < other.end:
            return True
        return False

    def __gt__(self, other):
        if self.start > other.start:
            return True
        if self.start < other.start:
            return False
        if self.end > other.end:
            return True
        return False

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end


class Interval(BaseInterval):
    def __init__(self, start, end):
        self._start = start
        self._end = end

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    def __repr__(self):
        return 'Interval({}, {})'.format(self.start, self.end)


class IntervalNode(object):
    def __init__(self, data, left=None, right=None):
        self.left = left
        self.right = right
        self.data = data

    def find_min(self):
        current_node = self
        while current_node.left_child:
            current_node = current_node.left_child
        return current_node

    def lookup(self, data, parent=None):
        if data < self.data:
            if self.left is None:
                return None, None
            return self.left.lookup(data, self)
        elif data > self.data:
            if self.right is None:
                return None, None
            return self.right.lookup(data, self)
        else:
            return self, parent

    def delete(self, data):
        # get node containing data
        node, parent = self.lookup(data)
        if node is not None:
            if node.is_leaf:
                # if node has no children, just remove it
                # check if it is not the root node
                if parent.left is node:
                    parent.left = None
                else:
                    parent.right = None
                del node
            elif node.left and node.right:
                # if node has 2 children
                # find its successor
                parent = node
                successor = node.right
                while successor.left:
                    parent = successor
                    successor = successor.left
                # replace node data by its successor data
                node.data = successor.data
                # fix successor's parent node child
                if parent.left == successor:
                    parent.left = successor.right
                else:
                    parent.right = successor.right
            else:
                # if node has 1 child
                # replace node by its child
                if node.left:
                    n = node.left
                else:
                    n = node.right

                if parent.left is node:
                    parent.left = n
                else:
                    parent.right = n
                del node
        else:
            raise LookupError()

    @property
    def is_leaf(self):
        return self.left is None and self.right is None

    def __iter__(self):
        for i in node_traverse(self):
            yield i