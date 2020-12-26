# Where libMesh is installed.
prefix              := /usr/local

# the libtool to use in these defined rules.
libmesh_LIBTOOL ?= ${prefix}/contrib/bin/libtool

libmesh_CXX      := $(shell ${prefix}/bin/libmesh-config --cxx)
libmesh_INCLUDE  := $(shell ${prefix}/bin/libmesh-config --include)
libmesh_CPPFLAGS := $(shell ${prefix}/bin/libmesh-config --cppflags)
libmesh_CXXFLAGS := $(shell ${prefix}/bin/libmesh-config --cxxflags)
libmesh_LIBS     := $(shell ${prefix}/bin/libmesh-config --libs)

all : mesh

mesh.o : mesh.cpp
	$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c $^ -o $@

mesh : convmesh.o mesh.o
	@$(libmesh_LIBTOOL) --tag=CXX --mode=link $(libmesh_CXX) $(libmesh_CXXFLAGS) -rdynamic $^ -o $@ $(libmesh_LIBS) -lmpfr -lgmpxx -lgmp /usr/lib/x86_64-linux-gnu/libCGAL.so.13.0.2 -lmpfr -lgmpxx -lgmp

clean :
	rm mesh.o mesh
