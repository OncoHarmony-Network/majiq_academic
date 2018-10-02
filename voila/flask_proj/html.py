class Html:
    def __init__(self):
        self._tags = []

    def __str__(self):
        return ''.join(self._str_helper())

    def _str_helper(self, idx=0):

        try:
            arr = self._tags[idx]
        except IndexError:
            return ''

        if len(arr) == 1:
            yield arr[0]
        else:
            yield arr[0]
            yield from self._str_helper(idx + 1)
            yield arr[1]

    def tag(self, tag_v, **kwargs):
        classes = kwargs.pop('classes', [])
        if classes:
            kwargs['class'] = ' '.join(classes)

        attrs = [k.replace('_', '-') + '="' + str(v) + '"' for k, v in kwargs.items()]
        attrs = ' '.join(attrs)
        self._tags.append(['<' + tag_v + ' ' + attrs + '>', '</' + tag_v + '>'])

    def text(self, txt):
        self._tags.append([txt])

    def reset(self):
        self._tags = []
