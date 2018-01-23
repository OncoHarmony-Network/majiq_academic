import multiprocessing


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