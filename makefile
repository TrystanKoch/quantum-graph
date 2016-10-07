########################################################################
#  quantum-graph Project Makefile 
#    - Compiles source for quantum graph modeling programs.
#    - Documents workflow for quantum graph modeling and data creation.
#    - Documents workflow for plot creation from data.
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
## quantum-graph
##
##  My intention with this project is to create efficient and useful
##  programs that solve eigenvalue problems on quantum graphs with a
##  variety of boundary conditions.
##  
##  Quantum graphs are collections of vertices, connected by bonds each
##  having a metric and a differential operator. Because of the
##  similarities of the Schrodinger equation, wave equation, and AC
##  circuit networks, solving the Helmholtz equation on the graph can
##  lead to solutions in each of these domains.
##
##  This project (will eventually) contain a common eigenvalue and 
##  eigenvector solver for the quantum graph problem, code allowing
##  the user to create various quantum graphs and boundary conditions
##  with relative ease, and code translating the base solutions into
##  physically meaningful eigenvalues and eigenfunctions for the 
##  quantum, wave, and electrical problems.
##
##  While I intend this mostly for my own use, once I finish using it
##  for my purposes, I wish to release it to others in the quantum graph
##  research community. Therefore, I have placed all of my code under
##  GPLv3, and used openly available free-software libraries.
##
##  In addition to using the makefile to compile my code, I have 
##  included a section containing my workflow. This section contains
##  many of the commands I used to create datasets, plots, and other
##  results. They are not intended to be perfect, or to be complete, but
##  the commands below will hopefully make it easier for others (and my
##  future self, honestly) to follow the steps I made and recreate my
##  work and results.
##
##  I originally created this code for my Ph.D. research at the
##  University of Maryland, College Park. It would not have been
##  possible without support from the Department of Physics, and The 
##  Institute for Research in Elecronics and Applied Physics. 
##
##  Additionally, Thomas Antonsen contributed much to this project as
##  my advisor at UMD.



########################################################################
### Directory variables

# Source file directory
SRCDIR=$(CURDIR)/src

# Header directory
INCDIR=$(CURDIR)/headers

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
### Definitions, targets, and paths

TARGETS=bound_qg_roots \
	find_bounded_qg_roots \
	take_differences \
	optimize_histogram


LIBRARIES=quantumgraph \
		quantumgraphobject \
		quantumgraphrootfinding


DIRECTORIES=$(OUTDIR) $(BUILDDIR) data


vpath %     $(OUTDIR)
vpath %.h   $(INCDIR)
vpath %.cpp $(SRCDIR)
vpath %.o   $(BUILDDIR)


CODEOBJECTS=$(TARGETS:%=%.o)
LIBRARYOBJECTS=$(LIBRARIES:%=%.o)



########################################################################
### Main recipe

# "make" creates all unpackaged directories and all executables  
code: | $(DIRECTORIES)
code: $(TARGETS) $(LIBRARYOBJECTS)



########################################################################
### Source File Compilation

.SECONDEXPANSION:

# Create Programs
$(TARGETS): $$(patsubst %,%.o,$$@) $(LIBRARYOBJECTS)
	@$(CC) $(CFLAGS) $(INCFLAGS) $(BUILDDIR)/$@.o \
		$(addprefix $(BUILDDIR)/,$(LIBRARYOBJECTS)) \
		-o $(OUTDIR)/$@ $(LDFLAGS) $(INCFLAGS)
	@echo "  Compiled $@"


# Make Object Files
$(CODEOBJECTS): $$(patsubst %.o,%.cpp,$$@)
	@$(CC) $(CFLAGS) $(INCFLAGS) \
		-c $(SRCDIR)/$(patsubst %.o,%.cpp,$@) \
		-o $(BUILDDIR)/$@ $(INCFLAGS)
	@echo "  Compiled $@"


# Create Class and Function Libraries
$(LIBRARYOBJECTS): $$(patsubst %.o,%.cpp,$$@) $$(patsubst %.o,%.h,$$@)
	@$(CC) $(CFLAGS) $(INCFLAGS) \
		-c $(SRCDIR)/$(patsubst %.o,%.cpp,$@) \
		-o $(BUILDDIR)/$@ 
	@echo "  Compiled $@"



########################################################################
### Directories

# For the directories that must be made by make.
$(DIRECTORIES):
	@mkdir -p $@
	@echo "  Made directory: $@/"
	

# Creates an archive file.
data/%:
	@mkdir -p $@
	@echo "  Made directory: $@/"



########################################################################
### Housekeeping

.PHONY: tidy clean cleandata veryclean


# "make tidy" removes only the temporary files created by editing
tidy: 
	@rm -rf *~ */*~ */*/*~ */*/*/*~
	@echo "  Removed temporary files"


# "make clean" removes the object files and any temp files.
clean: tidy
	@rm -rf $(BUILDDIR)
	@echo "  Removed object files"


# "make cleandata" removes data files to an archive directory.
# Names of files are determined by the date of cleaning.
# Files with the same name are backed up and not overwritten. 
# The archive files are 
cleandata: | data/$(shell date -I)
	@if mv --backup=t *.dat data/$(shell date -I) 2>/dev/null; \
		then echo "  Archived data to: data/`date -I`/"; fi


# "make veryclean" also removes the compiled target binaries.
# It runs both clean (by extension tidy) and cleandata, so all that is
# left in the folder are the source files.
veryclean: clean cleandata
	@rm -rf $(OUTDIR)
	@echo "  Removed program files"


