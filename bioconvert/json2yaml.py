"""Convert :term:`JSON` to :term:`YAML` format"""
import yaml, json
from bioconvert import ConvBase

__all__ = ["JSON2YAML"]


class JSON2YAML(ConvBase):
    """Convert :term:`JSON` file into :term:`YAML` file

    Conversion is based on yaml and json standard Python modules
    Indentation is set to 4 by default and affects the sections (not the list).
    For example::

        fruits_list:
        - apple
        - orange
        section1:
            do: true
            misc: 1

    """
    input_ext = [".json"]
    output_ext = [".yaml"]
    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input JSON file
        :param str outfile: input YAML file. 
        """
        super(JSON2YAML, self).__init__(infile, outfile, *args, **kargs)

    def get_yaml(self):
        with open(self.infile, "r") as infile:
            data = json.load(infile)
        yamldata = yaml.dump(data, Dumper=yaml.dumper.Dumper,
                default_flow_style="", indent=4)
        return yamldata

    def __call__(self):
        with open(self.outfile, "w") as outfile:
            outfile.write(self.get_yaml())

