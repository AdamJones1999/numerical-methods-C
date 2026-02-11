'''
OS: Ubuntu Linux
Script which does the following:
1. Read little endian binary data output from 
numerical methods in c into python to plot.
2. Convert binary64 in memoryview to double.
'''
import os
import struct

datapath = "../numerical-methods-C/test1.data"
size_bytes = os.path.getsize(datapath)

DOUBLE_SIZE = 8

with open(datapath, 'rb') as f: # "read binary" mode 
	buffer = bytearray(size_bytes)
	view = memoryview(buffer)
	nbytes = f.readinto(view) # binary mode required for this to work
	print(f"size of data in bytes: {size_bytes}")
	print(f"hex dump of data: {buffer.hex()}")


# converting binary64 to floating point using structs
n_doub = nbytes // DOUBLE_SIZE
print(f"n_doub: {n_doub}")

st = struct.Struct(f'<{n_doub}d')
double_arr = bytearray(n_doub)
conv_view = st.unpack(view)
print(f"converted double values: {conv_view}")

f.close()