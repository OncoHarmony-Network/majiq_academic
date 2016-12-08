from multiprocessing import Process, Pool, Manager
from multiprocessing.queues import JoinableQueue

from voila.constants import PROCESS_COUNT


class ProducerConsumer(object):
    def __init__(self):
        self.queue = None
        self.manager_dict = None

    def __enter__(self):
        raise NotImplementedError()

    def run(self):
        self.close()

        self.manager_dict = Manager().dict()
        self.queue = JoinableQueue()

        producer_proc = Process(target=self._producer)
        producer_proc.daemon = True
        producer_proc.start()

        pool = Pool(PROCESS_COUNT, self._worker)

        producer_proc.join()
        self.queue.join()

        pool.close()
        self.queue.close()
        producer_proc.terminate()

        self.__enter__()

    def _worker(self):
        raise NotImplementedError()

    def _producer(self):
        raise NotImplementedError()

    def close(self):
        raise NotImplementedError()
