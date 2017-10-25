"""
Possible Conversion
====================

Plot directed graph of possible conversion


"""
#################################################
#
from bioconvert.core.graph import create_graph


#####################################################
# Get a data set (BAM file) for testing
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 250
create_graph("conversion.png")

#####################################################
#  PNG does not work on RTD with graphviz from conda
#  SVG cannot be shown 
from pylab import imshow, imread, xticks, yticks, gca
imshow(imread("conversion.png"), interpolation="nearest")
xticks([])
yticks([])
ax = gca()
ax.axis('off')
