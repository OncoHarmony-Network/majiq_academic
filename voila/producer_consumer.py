from multiprocessing import Process, Pool, Manager
from multiprocessing.queues import JoinableQueue

from voila.constants import PROCESS_COUNT


class ProducerConsumer(object):
    def __init__(self):
        """
        Class to help with the producer consumer pattern.
        """
        self.queue = None
        self.manager_dict = None
        self.manager = None

    def __enter__(self):
        raise NotImplementedError()

    def manager_shutdown(self):
        """
        When the data has been copied from the manager instance, it can be shutdown.
        :return: None
        """
        self.manager.shutdown()

    def run(self):
        # this process can't have the same hdf5 file open twice
        self.close()

        self.manager = Manager()
        self.manager_dict = self.manager.dict()
        self.queue = JoinableQueue()

        producer_proc = Process(target=self._producer)
        producer_proc.start()

        pool = Pool(PROCESS_COUNT, self._worker)

        producer_proc.join()
        self.queue.join()

        pool.close()
        self.queue.close()

        # when finished, open file back up for this process
        self.__enter__()

    def _worker(self):
        """
        Function that the worker processes run.
        :return: None
        """
        raise NotImplementedError()

    def _producer(self):
        """
        Producer process that fill the queue.
        :return: None
        """
        raise NotImplementedError()

    def close(self):
        raise NotImplementedError()
