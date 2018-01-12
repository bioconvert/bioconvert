"""
Available methods per converter
=====================================

Plot number of implemented methods per converter.


"""
#################################################
#
from bioconvert.core.registry import Registry

mapper = Registry()

# The available unique converters
converters = set([mapper[this] for this in mapper._ext_registry])

# the number of methods per converter
data = [len(x.available_methods) for x in converters]

print("Number of converters: {}".format(len(converters)))
print("Number of methods : {}".format(sum(data)))

#####################################################
from pylab import hist, clf, xlabel, grid

clf()
hist(data, range(17), ec="k", zorder=2)
xlabel("Number of methods")
grid(zorder=-1)


