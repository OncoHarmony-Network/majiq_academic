from multiprocessing.pool import Pool


class VoilaPool:
    class __VoilaPool:
        def __init__(self, processes):
            self.processes = int(processes)
            self.pool = Pool(processes=self.processes)

        def __str__(self):
            return repr(self)

        def apply_async(self, *args, **kwargs):
            return self.pool.apply_async(*args, **kwargs)

    instance = None

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
