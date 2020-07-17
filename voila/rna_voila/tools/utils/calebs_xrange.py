# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


def calebs_xrange(start, end, by):
    """
    Does an xrange, but can handle floats instead of integers.

    Warning: floats don't behave quite as you'd expect on computers.
     Be careful how you use this function.
    """
    start = float(start)
    end = float(end)
    by = float(by)
    stuff_after_dec = str(1.0).split(".")[1]
    # Figure out which number has the smallest digits, use that
    # to generate the multiplier (to prevent rounding error)
    for number in start, end, by:
        if len(str(number).split(".")[1]) > len(stuff_after_dec):
            stuff_after_dec = str(number).split(".")[1]
    multyplier = float(pow(10, len(stuff_after_dec)))
    iteratable = list()
    # I add 'int(By*Multyplier)' to the End here because I like things to be explicitly
    #  defined. This function starts and ends exactly how you tell it, rather than
    #  Python's stupid way. Blah. Sorry, Caleb was cranky when he wrote this.
    for i in range(int(start * multyplier),
                   int(end * multyplier) + int(by * multyplier),
                   int(by * multyplier)):
        iteratable.append(float(i) / multyplier)
    return iteratable
