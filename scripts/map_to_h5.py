#!/usr/bin/python
#
# (c) 2018 Gergely Katona <gergely.katona@gu.se>
#
import iotbx.ccp4_map
import sys
import h5py

map_object1 = iotbx.ccp4_map.map_reader(sys.argv[1])

with h5py.File(sys.argv[2], "w") as f:
    f.create_dataset('map', data=map_object1.data.as_numpy_array())
    f.close()
    
