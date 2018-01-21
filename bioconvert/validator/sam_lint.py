import colorlog
_log = colorlog.getLogger(__name__)
import re

__all__ = ['SAMLint']


class SAMLint(object):
    """SAM validator

    ::

        from bioconvert.validator.sam_lint import SAMLint
        from bioconvert import bioconvert_data
        filename = bioconvert_data("test_measles.sam")
        SAMLint(filename).validate()

    """
    def __init__(self, filename):
        self.filename = filename

    def validate(self):
        # read line by line. Skip all lines starting with @
        # Check that there are 11 fields (compulsary) and
        # that following optional fields follow the TAG:TYPE:VALUE format
        with open(self.filename, "r") as fh:
            for lineno, line in enumerate(fh.readlines()):
                if line.startswith("@"):
                    continue
                fields = line.split()
                if len(fields) < 11:
                    msg = "Lines must contain at least 11 fields."
                    msg += "Found {} on line {}"
                    raise ValueError(msg.format(len(fields), lineno))
                if len(fields)>11:
                    for field in fields[11:]:
                        if re.match("(\S+):(\S+):(\S+)", field) is None:
                            msg = "Optional fields must follow the tag:type:value "
                            msg += "format. Found {} on line {}"
                            raise ValueError(msg.format(len(fields), lineno))





