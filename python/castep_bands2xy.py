#!/usr/bin/env python
################################################################################
#
# castep_bands2xy.py
#
# Converts a castep .bands file to a .bands_xy file.
#
# Usage: castep_bands2xy.py FILE
#
################################################################################
#
# Copyright 2013 Kane O'Donnell
#
#     This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# 
# NOTES
#
# 1. Configured for diamond unit cell.
#
################################################################################

from numpy import savetxt
import argparse
import esc_lib as el

parser = argparse.ArgumentParser(description="Convert a .bands file to a .bands_xy file where the columns are individual bands.")
parser.add_argument('filename', help="The .bands file to convert.")
args = parser.parse_args()  

bands, props = el.castep_read_bands(args.filename)
savetxt(args.filename+"_xy", bands)

