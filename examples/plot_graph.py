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
create_graph("conversion", format="eps")

#####################################################
#  PNG does not work on RTD with graphviz from conda
#  SVG cannot be shown 
from pylab import imshow, imread, xticks, yticks, gca
imshow(imread("conversion.eps"))
xticks([])
yticks([])
ax = gca()
ax.axis('off')
