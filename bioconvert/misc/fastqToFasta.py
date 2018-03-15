# -*- coding: utf-8 -*-
#!/opt/local/bin/python

import sys
import mmap

with open(sys.argv[1], "r+") as inp:
    with open(sys.argv[2], "wb") as out:
        mapp = mmap.mmap(inp.fileno(), 0)
        line = mapp.readline()
        while line:
            out.write(b">")
            out.write(line[1:])
            out.write(mapp.readline())
            mapp.readline()
            mapp.readline()
            line = mapp.readline()
        mapp.close()