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

        if attrs:
            open_tag = '<' + tag_v + ' ' + attrs + '>'
        else:
            open_tag = '<' + tag_v + '>'

        close_tag = '</' + tag_v + '>'

        self._tags.append([open_tag, close_tag])
        return self

    def text(self, txt):
        self._tags.append([txt])
        return self

    def children(self, *args):
        self._tags.append([''.join(str(a) for a in args)])
        return self

    def reset(self):
        self._tags = []
