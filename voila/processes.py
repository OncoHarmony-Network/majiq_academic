import multiprocessing

from voila.utils.voila_log import voila_log


class VoilaPool:
    class __VoilaPool:
        def __init__(self, processes):
            self.processes = int(processes)
            self.pool = multiprocessing.Pool(processes=self.processes)

        def __str__(self):
            return repr(self)

        def apply_async(self, *args, **kwargs):
            return self.pool.apply_async(*args, **kwargs)

    instance = None

    def __str__(self):
        return repr(self.instance)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __init__(self, processes=None):
        if not VoilaPool.instance:
            VoilaPool.instance = VoilaPool.__VoilaPool(processes)

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def close(self):
        self.pool.close()
        self.pool.join()


class VoilaQueue:
    class __VoilaQueue:
        def __init__(self, nprocs):
            self.Mgr = multiprocessing.Manager()
            self.queue = self.Mgr.Queue(maxsize=nprocs * 2)

            voila_log().debug('Manager PID {}'.format(self.Mgr._process.ident))

        def __str__(self):
            return repr(self)

    instance = None

    def __str__(self):
        return repr(self.instance)

    def __init__(self, filling_fn=None, nprocs=None):
        self.filling_fn = filling_fn
        if not VoilaQueue.instance:
            VoilaQueue.instance = VoilaQueue.__VoilaQueue(nprocs)

    def __enter__(self):
        q = self.instance.queue
        e = self.instance.Mgr.Event()

        self.process = multiprocessing.Process(target=self.filling_fn, args=(q, e))
        self.process.start()

        return q, e

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.instance.queue.empty():
            raise Exception('Queue is not empty yet...')

    def close(self):
        self.instance.Mgr.shutdown()
