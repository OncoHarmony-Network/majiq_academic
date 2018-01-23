import multiprocessing
import os
import time

from voila.api import Voila


class VoilaPool:
    class __VoilaPool:
        def __init__(self, nproccesses):
            self.pool = multiprocessing.Pool(processes=nproccesses, maxtasksperchild=1)

        def __str__(self):
            return repr(self)

    instance = None

    def __init__(self, nproccesses=None):
        if not VoilaPool.instance and nproccesses is None:
            raise Exception('Need to set the number of processes')

        if not VoilaPool.instance:
            VoilaPool.instance = VoilaPool.__VoilaPool(nproccesses)

    def __getattr__(self, name):
        return getattr(self.instance, name)


def chunkify(lst, n):
    lst = tuple(lst)
    for i in range(n):
        yield lst[i::n]


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def f(xs):
    with Voila(os.path.expanduser('~/Downloads/HL-60.psi.voila'), 'r') as vv:
        return tuple(vv.get_voila_lsv(x, y).strand for (x, y) in xs)


if __name__ == "__main__":
    processes = 4

    # val = [random.getrandbits(999999999) for x in range(25)]

    pool = VoilaPool(processes).pool

    ts = time.time()

    with Voila(os.path.expanduser('~/Downloads/HL-60.psi.voila'), 'r') as v:
        lsvs = tuple(v.get_lsvs())

    multiple_results = [pool.apply_async(f, (x,)) for x in chunkify(lsvs, processes)]
    tuple(res.get() for res in multiple_results)

    print(time.time() - ts)
