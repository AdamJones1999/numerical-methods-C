'''
OS: Ubuntu Linux
File to read little endian binary data output from 
numerical methods in c into python to plot.
'''
import os

datapath = "../numerical-methods-C/test.data"
size_bytes = os.path.getsize(datapath)


with open(datapath, 'rb') as f: # "read binary" mode 
	buffer = bytearray(size_bytes)
	view = memoryview(buffer)
	nbytes = f.readinto(view) # binary mode required for this to work
	print(buffer.hex())
	print(size_bytes)
f.close()