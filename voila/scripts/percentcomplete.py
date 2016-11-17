from datetime import datetime, timedelta


class PercentComplete(object):
    def __init__(self, length, output=True):
        self.length = length
        self.count = float(0)
        self.start_time = datetime.now()
        self.output = output

    def __nonzero__(self):
        is_complete = self.length == self.count
        if is_complete:
            self.end()
        else:
            self.print_percent()
        return not is_complete

    def print_percent(self):
        self.count += 1

        if self.output:
            now = datetime.now()
            delta = (now - self.start_time).total_seconds()
            percent = round((self.count / self.length) * 100, 3)
            total = timedelta(seconds=(self.length * delta) / self.count)
            remaining = self.start_time + total - now
            rate = round(self.count / delta, 1)

            print '\r', '{0}%'.format(percent),
            print '--', ' -- '.join([
                'Estimated Remaining Time: {0}'.format(remaining),
                'Jobs/Second: {0}'.format(rate)
            ]),

    def range(self):
        return range(self.length)

    @staticmethod
    def end():
        print ''
