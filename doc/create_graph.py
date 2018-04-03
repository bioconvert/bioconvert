from bioconvert.core.graph import create_graph
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 250
create_graph("conversion.png", color_for_disabled_converter="black")
create_graph("conversion-local.png")
