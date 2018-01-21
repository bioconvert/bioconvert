import colorlog
_log = colorlog.getLogger(__name__)


__all__ = ['GFALint']


class GFALint(object):
    """

    see https://github.com/sjackman/gfalint/
    """
    def __init__(self, filename):
        self.filename = filename


    def validate(self):
        # read line by line. Checks 
        # - lines start with HT, VT or ED
        # - lines must be tab delimited
        # - lines VT must have 2 fields only
        with open(self.filename, "r") as fh:
            for i, line in enumerate(fh.readlines()):
                if line[0] in "#EFGHLOPSU":
                    pass
                elif len(line.strip()) == 0:
                    _log.warning("Found empty line on line %s" % line)
                else:
                    raise ValueError("Unknown starting field (%s) on line %s" % (line[0], i))

