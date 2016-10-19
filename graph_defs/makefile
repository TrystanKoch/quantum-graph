########################################################################
#  quantum-graph definitions makefile 
#    - Compiles source for quantum graph definition-outputting programs
#
#  Copyright (C) 2015-2016 Trystan Koch
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


########################################################################
## quantum-graph definitions
##
##  These programs contain definitions for quantum graphs that I have
##  created by hand for my research. Because the graphs can get 
##  complicated, these programs use QuantumGraphObject to create a graph
##  by defining nodes and bonds and connecting them. By running the code
##  we output a file containing the graph's scattering matrix and bond
##  lengths. This file may then be directed into certain other programs
##  like the bound_qg_roots program.




########################################################################
### Directory variables

# Source file directory
SRCDIR=$(CURDIR)/src

# Header directory
INCDIR=$(CURDIR)/../headers

# Place where the library object files are
LIBOBJDIR=$(CURDIR)/../obj

# Output object files to here
BUILDDIR=$(CURDIR)/obj

# Executable file destination
OUTDIR=$(CURDIR)/bin



########################################################################
### Compiler Flag definitions

# Uses the gcc C++ compiler.
CC=g++


# Uses the gcc C++ linker.
LD=g++


# -DLAPACK_COMPLEX_CUSTOM
#     Lapack complex type is redefined to the gsl_complex type in code.
# -std=c++11 
#     Use the features of C++11.
# -Wall
#     Show all warnings in the compilation. Always fix these.
# -O3
#     Optimize code using all possible means of optimization. 
CFLAGS=-DLAPACK_COMPLEX_CUSTOM -std=c++11  -Wall -O3


# -lgsl
#     Link the GNU Scientific Library.
# -lgslcblas
#     Link the GSL libraries for BLAS computation.
# -llapacke
#     Link the C version of LAPACK linear algebra package.
# -L$(INCDIR)
#     Link the headers from the headers folder
LDFLAGS=-lgsl -lgslcblas -llapacke -L$(INCDIR)


# -I/usr/include/gsl
#     Include the GNU Scientific Library.
# -I$(INCDIR)
#     Include the file headers.
INCFLAGS=-I/usr/include/gsl -I$(INCDIR)




########################################################################
### Code compilation
CODETARGETS=neuman_complete_five_graph \
		neuman_3bond_star \
		neuman_tetrahedron_input_loss \
		neuman_tetrahedron_input_loss_two_port \
		neuman_tetrahedron_nonreciprocal_lines \
		quotient_graph_attempt_1 \
		quotient_graph_2 \
		quotient_graph_3 \
		quotient_graph_5 \
		quotient_graph_6


LIBRARIES=quantumgraph \
	quantumgraphobject \


DIRECTORIES=$(OUTDIR) $(BUILDDIR) data


vpath %     $(OUTDIR)
vpath %.h   $(INCDIR)
vpath %.cpp $(SRCDIR)
vpath %.o   $(BUILDDIR)
vpath %.o   $(LIBOBJDIR)


CODEOBJECTS=$(CODETARGETS:%=%.o)
LIBRARYOBJECTS=$(LIBRARIES:%=%.o)


# "make" will create the program files only
code: $(CODETARGETS) 


########################################################################
### Source File Compilation

.SECONDEXPANSION:

# Create Programs
$(CODETARGETS): $$(patsubst %,%.o,$$@) $(LIBRARYOBJECTS)
	@$(CC) $(CFLAGS) $(INCFLAGS) $(BUILDDIR)/$@.o \
		$(addprefix $(LIBOBJDIR)/,$(LIBRARYOBJECTS)) \
		-o $(OUTDIR)/$@ $(LDFLAGS) $(INCFLAGS)
	@echo "  Compiled $@"
	
# Create Object Files
$(CODEOBJECTS): $$(patsubst %.o,%.cpp,$$@)
	@$(CC) $(CFLAGS) $(INCFLAGS) \
		-c $(SRCDIR)/$(patsubst %.o,%.cpp,$@) \
		-o $(BUILDDIR)/$@ $(INCFLAGS)
	@echo "  Compiled $@"
	



########################################################################
### Housekeeping

.PHONY: tidy clean cleandata veryclean

# "make tidy" removes only the temporary files created by editing
tidy: 
	@rm -f *~ */*~
	@echo "  Removed temporary files"


# "make clean" removes the object files and any temp files.
clean: tidy
	@rm -f $(BUILDDIR)/*.o
	@echo "  Removed object files"


# "make cleandata" removes data files to an archive directory.
# Files with the same name are backed up and not overwritten. 
# The archive files are 
cleandata: | data/$(shell date -I)
	@if test -e *.dat; \
		then mv --backup=t *.dat data/$(shell date -I) \
		&& echo "  Archived data files to: data/`date -I`/"; fi

# Creates an archive file with the present date as the name
data/%:
	@mkdir -p $@
	@echo "  Made directory: $@/"


# "make veryclean" also removes the compiled target binaries.
# It runs both clean (by extension tidy) and cleandata, so all that is
# left in the folder are the source files.
veryclean: clean
	@rm -f $(OUTDIR)/*
	@echo "  Removed program files"

