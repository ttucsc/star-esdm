#!/usr/bin/env python
#
#   python script to output the canonical path to a given file.
#
#   NOTE:  returns a path even if the file doesn't exist...just resolves all valid paths along the way.
#
import os
import sys

print(os.path.realpath(sys.argv[1]))
