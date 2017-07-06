from voila.tools import Tool


class DoesNotMatter(Tool):
    help = 'example code who cares'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('--yup', action='store_true', help='this is yup')
        parser.add_argument('--upper', action='store_true', help='upper case')
        return parser

    def run(self, args):
        print(hello_world(capitalize=args.upper))


def hello_world(capitalize=False):
    if capitalize:
        return 'Hello, world!'.upper()
    else:
        return 'Hello, world!'
