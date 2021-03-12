"""
majiq_args.py

Helper classes/functions for checking valid arguments directly at command-line
"""

import argparse
from typing import (
    Any,
    Callable,
    Optional,
    Set,
    Union,
)


def check_characters_factory(
    valid_characters: Set[str],
    valid_description: Optional[str] = None,
) -> Callable[[Any], str]:
    """Create function to raise argparse type error or return string

    Create function to raise argparse.ArgumentTypeError if input contains
    invalid characters, otherwise, return str(x)

    Parameters
    ----------
    valid_characters: Set[str]
        characters that are acceptable for this restricted string
    valid_description: Optional[str]
        description of acceptable characters in input
    """

    def check_characters(x: Any) -> str:
        x = str(x)
        set_x = set(x)
        if not set_x.issubset(valid_characters):
            invalid_chars = set_x.difference(valid_characters)
            prefix = (
                f"may only contain {valid_description}: " if valid_description else ""
            )
            raise argparse.ArgumentTypeError(
                f"{prefix}{x} contains invalid characters {invalid_chars}"
            )
        return x

    return check_characters


def check_nonnegative_factory(
    cast_fn: Callable[[Any], Union[int, float]], reject_zero: bool
):
    """Returns argparse type function to casted type, rejecting negative values"""
    invalid_predicate = (lambda x: x <= 0) if reject_zero else (lambda x: x < 0)
    invalid_message = "parameter must be {}, {{}} is not".format(
        "positive" if reject_zero else "nonnegative"
    )

    def return_function(value):
        ivalue = cast_fn(value)
        if invalid_predicate(ivalue):
            raise argparse.ArgumentTypeError(invalid_message.format(value))
        return ivalue

    return return_function


def StoreRequiredUniqueActionFactory():
    """Return class that pools shared values (do not allow overlap)"""

    class StoreRequiredUniqueAction(argparse.Action):

        # class variable with shared values known to all instances
        shared_values = {}  # Dict[str, Set[str]]

        @staticmethod
        def repeated_values(values):
            """Returns repeated values (that occur more than twice)

            Returns
            -------
            Set[str]
            """
            unique_values = set(values)  # set of unique values
            repeated_values = set()  # set of values that were repeated
            if isinstance(values, list) and len(values) != len(unique_values):
                # values are not all unique...
                repeated_values = set()  # determine which values were repeated
                # remove unique values (if done twice, added to repeated_values)
                for x in values:
                    try:
                        unique_values.remove(x)
                    except KeyError:
                        repeated_values.add(x)
            return repeated_values

        @classmethod
        def overlapping_values(cls, dest, values):
            """Updates shared_values and returns overlapping values with other
            parameters

            Parameters
            ----------
            dest: str
                Key to ignore in shared_values
            values:
                Values to compare with other shared values

            Returns
            -------
            Dict[str, Set[str]]
                Keys for other parameter names, values correspond to overlaps
                between them
            """
            unique_values = set(values)
            overlapping_values = {
                other_key: other_values & unique_values
                for other_key, other_values in cls.shared_values.items()
                if other_key != dest and other_values & unique_values
            }
            # update shared values
            cls.shared_values[dest] = unique_values
            # remove any keys with zero overlaps
            return {k: v for k, v in overlapping_values.items() if v}

        def __call__(self, parser, namespace, values, option_string=None):
            """Check that values are unique if list, not in shared_values"""
            # no repeated or overlapping values
            repeated_values = self.repeated_values(values)
            if repeated_values:
                raise argparse.ArgumentError(
                    self, f"values must be unique (repeated values: {repeated_values})"
                )
            overlapping_values = self.overlapping_values(self.dest, values)
            if overlapping_values:
                raise argparse.ArgumentError(
                    self,
                    "values must not overlap shared parameters"
                    f" (overlaps: {overlapping_values})",
                )
            # no non-unique or overlapping values, so save updated values
            setattr(namespace, self.dest, values)

    return StoreRequiredUniqueAction
