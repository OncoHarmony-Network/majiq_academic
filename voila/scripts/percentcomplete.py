from datetime import datetime, timedelta


class PercentComplete(object):
    def __init__(self, length, output=True):
        self.length = length
        self.count = float(0)
        self.start_time = datetime.now()
        self.output = output

    def print_percent(self):
        self.count += 1

        now = datetime.now()
        delta = (now - self.start_time).total_seconds()
        percent = round((self.count / self.length) * 100, 3)
        total = timedelta(seconds=(self.length * delta) / self.count)
        remaining = self.start_time + total - now
        rate = round(self.count / delta, 1)

        if self.output:
            print '\r', '{0}%'.format(percent),
            if self.count > 10:
                print '--', ' -- '.join([
                    'Estimated Remaining Time: {0}'.format(remaining),
                    'Jobs/Second: {0}'.format(rate)
                ]),

    def range(self):
        return range(self.length)

    @staticmethod
    def end():
        print ''
