########################################################################
#  quantum-graph definitions makefile 
#    - Compiles source for quantum graph definition-outputting programs
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
## quantum-graph definitions
##
##  Because various 



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
CODETARGETS=create_neumann_complete_five_graph \
		create_neumann_3bond_star \
		create_neumann_tetrahedron_input_loss \
		create_neumann_tetrahedron_input_loss_two_port

QGOBJECT=../quantumgraphobject.o

# "make" will create the program files only
code: $(CODETARGETS) 


########################################################################
### Source File Compilation

### Create Programs

create_neumann_tetrahedron_input: create_neumann_tetrahedron_input.o \
				../quantumgraph.o $(QGOBJECT)
	@$(CC) $(CFLAGS) create_neumann_tetrahedron_input.o \
		../quantumgraph.o -o create_neumann_tetrahedron_input \
		$(LDFLAGS)
	@echo "  Compiled $@"

create_neumann_tetrahedron_input_loss: create_neumann_tetrahedron_input_loss.o \
				../quantumgraph.o $(QGOBJECT)
	@$(CC) $(CFLAGS) create_neumann_tetrahedron_input_loss.o \
		../quantumgraph.o $(QGOBJECT) -o create_neumann_tetrahedron_input_loss \
		$(LDFLAGS)
	@echo "  Compiled $@"
	
create_neumann_tetrahedron_input_loss_two_port: create_neumann_tetrahedron_input_loss_two_port.o \
				../quantumgraph.o $(QGOBJECT)
	@$(CC) $(CFLAGS) create_neumann_tetrahedron_input_loss_two_port.o \
		../quantumgraph.o $(QGOBJECT) -o create_neumann_tetrahedron_input_loss_two_port \
		$(LDFLAGS)
	@echo "  Compiled $@"

create_neumann_complete_five_graph: create_neumann_complete_five_graph.o \
				../quantumgraph.o $(QGOBJECT)
	@$(CC) $(CFLAGS) create_neumann_complete_five_graph.o \
		../quantumgraph.o $(QGOBJECT) -o create_neumann_complete_five_graph \
		$(LDFLAGS)
	@echo "  Compiled $@"

create_neumann_3bond_star: create_neumann_3bond_star.o \
				../quantumgraph.o $(QGOBJECT)
	@$(CC) $(CFLAGS) create_neumann_3bond_star.o \
		../quantumgraph.o $(QGOBJECT) -o create_neumann_3bond_star \
		$(LDFLAGS)
	@echo "  Compiled $@"


### Make Object Files

create_neumann_tetrahedron_input.o: create_neumann_tetrahedron_input.cpp
	@$(CC) $(CFLAGS) -c create_neumann_tetrahedron_input.cpp
	@echo "  Compiled $@"

create_neumann_tetrahedron_input_loss.o: create_neumann_tetrahedron_input_loss.cpp
	@$(CC) $(CFLAGS) -c create_neumann_tetrahedron_input_loss.cpp
	@echo "  Compiled $@"

create_neumann_tetrahedron_input_loss_two_port.o: create_neumann_tetrahedron_input_loss_two_port.cpp
	@$(CC) $(CFLAGS) -c create_neumann_tetrahedron_input_loss_two_port.cpp
	@echo "  Compiled $@"
	
create_neumann_complete_five_graph.o: create_neumann_complete_five_graph.cpp
	@$(CC) $(CFLAGS) -c create_neumann_complete_five_graph.cpp
	@echo "  Compiled $@"

create_neumann_3bond_star.o: create_neumann_3bond_star.cpp
	@$(CC) $(CFLAGS) -c create_neumann_3bond_star.cpp
	@echo "  Compiled $@"


### Create Class and Function Libraries

../quantumgraph.o: ../quantumgraph.cpp ../quantumgraph.h
	@$(CC) $(CFLAGS) -c ../quantumgraph.cpp -o ../quantumgraph.o
	@echo "  Compiled $@"

../quantumgraphobject.o: ../quantumgraphobject.cpp ../quantumgraphobject.h
	@$(CC) $(CFLAGS) -c ../quantumgraphobject.cpp -o ../quantumgraphobject.o
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
veryclean: clean
	@rm -f $(CODETARGETS)
	@echo "  Removed program files"

