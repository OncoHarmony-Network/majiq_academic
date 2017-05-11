import signal
from multiprocessing import Process, Pool, Manager
from multiprocessing.queues import JoinableQueue

from voila import constants
from voila.utils.voila_log import voila_log


def voila_exit(error_no=1):
    voila_log().warning("Voila exiting...")
    exit(error_no)


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
        """
        Runs multiple process producer/consumer pattern.
        :return: None
        """

        # this process can't have the same hdf5 file open twice
        self.close()

        self.manager = Manager()
        self.manager_dict = self.manager.dict()
        self.queue = JoinableQueue()

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

        producer_proc = Process(target=self._producer)
        pool = Pool(constants.PROCESS_COUNT, self._worker)

        signal.signal(signal.SIGINT, original_sigint_handler)

        try:

            producer_proc.start()

            producer_proc.join()
            self.queue.join()

            pool.close()
            self.queue.close()

        except (KeyboardInterrupt, SystemExit):
            producer_proc.terminate()
            pool.terminate()
            self.manager_shutdown()
            voila_exit()

        # when finished, open file back up for this process
        self.__enter__()

    def dict(self, key, value):
        """
        Manager dictionary to hold data created by worker processes.

        Note: EOF and IO Error will happen when script is terminated early.
        :param key: key
        :param value: value
        :return: None
        """
        try:
            self.manager_dict[key] = value
        except EOFError:
            pass
        except IOError:
            pass

    def get_values(self):
        try:
            return self.manager_dict.values()
        except KeyboardInterrupt:
            voila_exit()

    def get_dict(self):
        try:
            return dict(self.manager_dict)
        except KeyboardInterrupt:
            voila_exit()

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
