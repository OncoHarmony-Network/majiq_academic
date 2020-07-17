from rna_voila.tools import Tool
from rna_voila.exceptions import VoilaException


class ExampleException(VoilaException):
    """

    Using VoilaException for your exception classes will allow them to work more seamlessly with the debug argument.

    """

    pass


class CtrlPsiDictionaryTool(Tool):
    """

    Add this class to run your program.  The name of the class doesn't actually matter, but it helps when it's
    descriptive.

    """

    help = 'This is an example of how to integrate your script into Voila Tools.'

    def arguments(self):
        """

        Collects the arguments required to run the program.

        :return: parser object
        """
        parser = self.get_parser()
        parser.add_argument('--upper', action='store_true', help='Print text in upper case')
        parser.add_argument('--raise-exception', action='store_true', help='Raises Voila exception')
        return parser

    def run(self, args):
        """

        Run the program using the command line arguments.

        :param args: command line arguments
        :return: None
        """

        if args.raise_exception:
            raise ExampleException('Forced Voila exception')

        example_text = hello_world(upper=args.upper)
        print(example_text)


def hello_world(upper=False):
    """

    For the purpose of this example, this is our program.

    :param upper: return text with only uppercase letters.
    :return: String
    """
    example_text = 'Hello, world!'
    if upper:
        return example_text.upper()
    return example_text
