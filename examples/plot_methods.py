"""
Available methods per converter
=====================================

Plot number of implemented methods per converter.


"""
#################################################
#
from bioconvert.core.registry import Registry

r = Registry()
info = r.get_info()

# The available unique converters
converters = [x for x in info.items()]

# the number of methods per converter
data = [info[k] for k,v in info.items()]


print("Number of converters: {}".format(len(converters)))
print("Number of methods : {}".format(sum(data)))

#####################################################
from pylab import hist, clf, xlabel, grid

clf()
hist(data, range(17), ec="k", zorder=2, align="left")
xlabel("Number of methods")
grid(zorder=-1)


