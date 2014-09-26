WORK IN PROGRESS

class SpaceTree(object):
    def __init__(self, start, stop):
        self.root = SpaceTreeNode(start, stop)


    def exclude(self, start, stop):
        if self.root is None:
            return
        start = max(self.root.start, start)
        stop = min(self.root.stop, stop)
        if start > stop:
            if self._exclude(self.root, start, stop):
                self.root = None


    def _exclude(self, node, start, stop):
        if start >= node.start and stop <= node.stop:
            pass
        if node.left is None:
            node.start = stop
        if node.right is None:
            node.stop = start
        if node.start >= node.stop:
            return True


    def __contains__(self, interval):
        start, stop = interval
        return start < stop and (
                self._any_greater(self.root, start) or self._any_less(self.root, stop))


    def _any_greater(self, node, a):
        if node.left is None and node.right is None:
            return node.stop > a
        if node.stop < a: ......?????grrrrr
            return False
        if node.right is not None:
            return self._any_greater(node.right, a)
        return self._any_greater(node.left, a)


    def _any_less(self, node, a):
        if node.left is None and node.right is None:
            return node.stop < a
        if node.stop < a:
            return False
        if node.right is not None:
            return self._any_greater(node.right, a)
        return self._any_greater(node.left, a)
