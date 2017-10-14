"""Convert :term:`YAML` to :term:`JSON` format"""
import yaml, json
from bioconvert import ConvBase

__all__ = ["YAML2JSON"]


class YAML2JSON(ConvBase):
    """Convert :term:`YAML` file into :term:`JSON` file

    Conversion is based on yaml and json standard Python modules

    .. note:: YAML comments will be lost in JSON output


    :reference: http://yaml.org/spec/1.2/spec.html#id2759572
    """
    input_ext = [".yaml"]
    output_ext = [".json"]

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input YAML file. 
        :param str outfile: input JSON file
        """
        super(YAML2JSON, self).__init__(infile, outfile, *args, **kargs)

    def get_json(self):
        """Return the JSON dictionary corresponding to the YAML input"""
        data = yaml.load(open(self.infile, "r"))
        return json.dumps(data, sort_keys=True, indent=4)

    def __call__(self):
        with open(self.outfile, "w") as outfile:
            outfile.write(self.get_json())

