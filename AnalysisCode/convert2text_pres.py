#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script simply takes care of parsing netcdf TC files and turning them in
to text files. Note that in this version we write out basin, too
"""
import TC_Utils as tc
import GeneralFunctions as GF
from netCDF4 import Dataset

ncfile="/media/gytm3/WD12TB/TropicalCyclones/TC-DeadlyHeat/Data/IBTrACS.ALL.hotel1.nc"
tco="/media/gytm3/WD12TB/TropicalCyclones/TC-DeadlyHeat/Data/Hurricanes_press.txt"
tc.nc2text_pres(ncfile=ncfile,thresh=980,fo=tco)



