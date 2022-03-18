# create the sphinx code to include all converters
# author: Thomas Cokelaer 2018
from bioconvert.core.registry import Registry

reg = Registry()


convs = list(reg.get_info())
names = sorted([l.__module__ for l in convs])


def underline(text, character="-"):
    return text + "\n" + character * len(text) + "\n"


print(underline("Reference converters", "="))
print(underline("Summary", "-"))
print(".. autosummary::\n")
for name in names:
    print("\t{}".format(name))


print(underline("All converters documentation", "-"))
print()
for name in names:
    print(
        """
.. automodule:: {}
    :members:
    :synopsis:
    :private-members:""".format(
            name
        )
    )
