# Where libMesh is installed.
prefix              := /usr/local

# the libtool to use in these defined rules.
libmesh_LIBTOOL ?= ${prefix}/contrib/bin/libtool

libmesh_CXX      := $(shell ${prefix}/bin/libmesh-config --cxx)
libmesh_INCLUDE  := $(shell ${prefix}/bin/libmesh-config --include)
libmesh_CPPFLAGS := $(shell ${prefix}/bin/libmesh-config --cppflags)
libmesh_CXXFLAGS := $(shell ${prefix}/bin/libmesh-config --cxxflags)
libmesh_LIBS     := $(shell ${prefix}/bin/libmesh-config --libs)

all : introduction_ex4

exact_solution.o : exact_solution.C
	$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c $^ -o $@

introduction_ex4.o : introduction_ex4.C
	$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c $^ -o $@

introduction_ex4 : introduction_ex4.o exact_solution.o
	@$(libmesh_LIBTOOL) --tag=CXX --mode=link $(libmesh_CXX) $(libmesh_CXXFLAGS) $^ -o $@ $(libmesh_LIBS)

clean :
	rm introduction_ex4.o introduction_ex4
