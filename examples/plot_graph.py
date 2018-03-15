"""
Possible Conversion
====================

Plot directed graph of possible conversion


"""
#################################################
#
from bioconvert.core.graph import create_graph



#####################################################
# If you use pygraphviz, you can have a good quality 
# image using:
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 250

############################################################
# In order to create the following image, you need graphviz
# and pygraphviz.
# If you cannot install those packages, you may use a singularity
# image like in the following example by setting the use_singularity
# parameter to True. This would work under Linux. Not tested on other systems
# yet.
create_graph("conversion.png", use_singularity=True)

#####################################################
#
from pylab import imshow, imread, xticks, yticks, gca
imshow(imread("conversion.png"), interpolation="nearest")
xticks([])
yticks([])
ax = gca()
ax.axis('off')
