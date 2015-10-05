########################################################################
#  quantum-graph Project Makefile 
#    - Compiles source for quantum graph modeling programs.
#    - Documents workflow for quantum graph modeling and data creation.
#    - Documents workflow for plot creation from data.
#
#  Copyright (C) 2015 Trystan Koch
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
### Compiler Flag definitions

# Uses the gcc C++ compiler.
CC=g++

# -DLAPACK_COMPLEX_CUSTOM
#     Lapack complex type is redefined to the gsl_complex type in code.
# -std=c++11 
#     Use the features of C++11.
# -I/usr/include/gsl
#     Include the GNU Scientific Library.
# -Wall
#     Show all warnings in the compilation. Always fix these.
# -O3
#     Optimize code using all possible means of optimization. 
CFLAGS=-DLAPACK_COMPLEX_CUSTOM -std=c++11 -I/usr/include/gsl -Wall -O3

# -lgsl
#     Link the GNU Scientific Library.
# -lgslcblas
#     Link the GSL libraries for BLAS computation.
# -llapacke
#     Link the C version of LAPACK linear algebra package.
LDFLAGS=-lgsl -lgslcblas -llapacke


########################################################################
### Code compilation
CODETARGETS=bound_qg_roots find_bounded_qg_roots 

# "make" will create the program files only
code: $(CODETARGETS) quantumgraphobject.o


########################################################################
### Workflow
RESULTTARGETS=roots.dat

#
results: $(RESULTTARGETS)


### Use these directives to calculate the roots of the quantum graph
### which is defined in bound_qg_roots.cpp (to change)

roots.dat: bound_qg_roots find_bounded_qg_roots
	time ./bound_qg_roots | ./find_bounded_qg_roots >roots.dat


########################################################################
### Source File Compilation

### Create Programs

bound_qg_roots: bound_qg_roots.o quantumgraph.o \
			quantumgraphrootfinding.o
	@$(CC) $(CFLAGS) bound_qg_roots.o quantumgraph.o \
		quantumgraphrootfinding.o -o bound_qg_roots $(LDFLAGS)
	@echo "  Compiled $@"

find_bounded_qg_roots: find_bounded_qg_roots.o quantumgraph.o \
			quantumgraphrootfinding.o
	@$(CC) $(CFLAGS) find_bounded_qg_roots.o quantumgraph.o \
		quantumgraphrootfinding.o -o find_bounded_qg_roots \
		$(LDFLAGS)
	@echo "  Compiled $@"




### Make Object Files

bound_qg_roots.o: bound_qg_roots.cpp
	@$(CC) $(CFLAGS) -c bound_qg_roots.cpp
	@echo "  Compiled $@"

find_bounded_qg_roots.o: find_bounded_qg_roots.cpp
	@$(CC) $(CFLAGS) -c find_bounded_qg_roots.cpp
	@echo "  Compiled $@"




### Create Class and Function Libraries

quantumgraph.o: quantumgraph.cpp quantumgraph.h
	@$(CC) $(CFLAGS) -c quantumgraph.cpp
	@echo "  Compiled $@"

quantumgraphobject.o: quantumgraphobject.cpp quantumgraphobject.h
	@$(CC) $(CFLAGS) -c quantumgraphobject.cpp
	@echo "  Compiled $@"

quantumgraphrootfinding.o: quantumgraphrootfinding.cpp \
			quantumgraphrootfinding.h
	@$(CC) $(CFLAGS) -c quantumgraphrootfinding.cpp
	@echo "  Compiled $@"


########################################################################
### Housekeeping

.PHONY: tidy clean cleandata veryclean

# "make tidy" removes only the temporary files created by editing
tidy: 
	@rm -f *~
	@echo "  Removed temporary files"


# "make clean" removes the object files and any temp files.
clean: tidy
	@rm -f *.o
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
veryclean: clean cleandata
	@rm -f $(CODETARGETS)
	@echo "  Removed program files"

