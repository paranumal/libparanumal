
#0: h_foo host pointer
#1: occa::memory o_foo = device.malloc(Nbytes, h_foo);
#2: o_foo.copyFrom(h_foo);
#3: o_foo.copyFrom(h_foo, Nbytes, offsetBytes);
#3: o_foo.copyFrom(h_foo, Nbytes, offsetBytes, "async: true"); // requires h_foo is pinned
#4: o_foo.copyTo(h_foo); 

