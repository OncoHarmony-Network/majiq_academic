import tblib.pickling_support
tblib.pickling_support.install()

from multiprocessing import Pool
import sys


class ExceptionWrapper(object):

    def __init__(self, ee):
        self.ee = ee
        __,  __, self.tb = sys.exc_info()

    def re_raise(self):
        raise self.ee.with_traceback(self.tb)
        # for Python 2 replace the previous line by:
        # raise self.ee, None, self.tb


# example how to use ExceptionWrapper

def inverse(i):
    """will fail for i == 0"""
    try:
        return 1.0 / i
    except Exception as e:
        return ExceptionWrapper(e)


def main():
    p = Pool(1)
    results = p.map(inverse, [1, 2, 0, 3])
    for result in results:
        if isinstance(result, ExceptionWrapper):
            result.re_raise()


if __name__ == "__main__":
    main()