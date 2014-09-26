from linked_list import LinkedList
import numpy as np

# XXX: LinkedList represents the list and the nodes at the same time... weird
# -> cannot add pointer to the last node (for fast append) as it would add the pointer (memory) to all nodes
# ...but does not use the self.data (head data) at all!? (because of remove?)

class IntervalLinkedList(LinkedList):
    """
    Store non-overlapping intervals in a linked list.

    The list is sorted.
    """

    def __init__(self, iterable=(), data=None, next=None):
        stop = np.nan
        for i in iterable:
            assert len(i) >= 2
            assert i[0] < i[1]
            assert not i[0] < stop
            stop = i[1]
        super(IntervalLinkedList, self).__init__(iterable=iterable, data=data, next=next)


    def __contains__(self, data):
        """
        Check whether the interval data is completely inside one of the intervals.
        """
        current = self.next
        while current:
            if (data[0] >= current.data[0] and data[0] < current.data[1]) and \
                    (data[1] > current.data[0] and data[1] <= current.data[1]):
                return True
            current = current.next
        return False


    def remove(self, data):
        """
        Remove the interval data from the list.

        The interval does not have to explicitly appear in the list.
        """
        current = self.next
        start = self
        # 1. deal with the start of the removed interval
        while current is not None:
            if data[0] < current.data[1]:
                if data[0] < current.data[0]: # data start is before current
                    break
                if data[0] >= current.data[0]: # data start is in current
                    if data[1] < current.data[1]: # data is totally in current
                        if data[0] > current.data[0]: # break current in two
                            new = LinkedList(None, (data[1], current.data[1]), current.next)
                            current.next = new
                            current.data = (current.data[0], data[0])
                        else: # data start is at current start (just shorten current)
                            current.data = (data[1], current.data[1])
                        return
                    start = current
                    current.data = (current.data[0], data[0])
                break
            start = current
            current = current.next
        if current is None:
            return
        # 2. deal with the end of the removed interval
        while current is not None:
            if data[1] < current.data[1]: # current is longer than data
                break
            current = current.next
        if current is None:
            start.next = None
            return
        if current != start:
            start.next = current
        if data[1] > current.data[0]:
            current.data = (data[1], current.data[1])


def test_in():
    x = IntervalLinkedList([(1, 10), (20, 30)])
    assert (0, 2) not in x
    assert (10, 21) not in x
    assert (9, 10) in x
    assert (5, 25) not in x
    assert (22, 30) in x
    assert (0, 1) not in x
    assert (30, 40) not in x
    assert (1, 10) in x
    assert (20, 30) in x
    assert (20, 31) not in x


def test_remove():
    x = IntervalLinkedList([(1, 1000)])
    assert str(x) == '(1, 1000)', str(x)
    x.remove((-10, 1))
    assert str(x) == '(1, 1000)', str(x)
    x.remove((1000, 2000))
    assert str(x) == '(1, 1000)', str(x)
    x.remove((20, 50))
    assert str(x) == '(1, 20)->(50, 1000)', str(x)
    x.remove((15, 25))
    assert str(x) == '(1, 15)->(50, 1000)', str(x)
    x.remove((40, 60))
    assert str(x) == '(1, 15)->(60, 1000)', str(x)
    x.remove((10, 70))
    assert str(x) == '(1, 10)->(70, 1000)', str(x)
    x.remove((10, 11))
    assert str(x) == '(1, 10)->(70, 1000)', str(x)
    x.remove((70, 71))
    assert str(x) == '(1, 10)->(71, 1000)', str(x)
    x.remove((80, 81))
    assert str(x) == '(1, 10)->(71, 80)->(81, 1000)', str(x)
    x.remove((90, 100))
    assert str(x) == '(1, 10)->(71, 80)->(81, 90)->(100, 1000)', str(x)
    x.remove((150, 180))
    assert str(x) == '(1, 10)->(71, 80)->(81, 90)->(100, 150)->(180, 1000)', str(x)
    x.remove((95, 190))
    assert str(x) == '(1, 10)->(71, 80)->(81, 90)->(190, 1000)', str(x)
    for i in range(200, 800, 10):
        x.remove((i, i+2))
    x.remove((200, 900))
    assert str(x) == '(1, 10)->(71, 80)->(81, 90)->(190, 200)->(900, 1000)', str(x)
    x.remove((0, 1001))
    assert str(x) == ''
