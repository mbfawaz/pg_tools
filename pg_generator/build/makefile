# Project name
PRJ_NAME = pg_gen

# Include paths other than the default
INCPATHS = -I/usr/include/suitesparse -I../include

# Libraries paths other than the default
LIBPATHS = -L/usr/lib/libblas.so

# Link libraries
LIBS = -lm -lumfpack -lboost_timer-mt -lboost_system -lboost_thread-mt -lamd -lcolamd -lcholmod -lccolamd -lcamd -lbtf -lcsparse -lcxsparse -lklu -lldl -lspqr -llapack -lblas 
LDOPT = -Wl,-rpath-link, -Wl,-rpath,
  
#
# Specify compiler
CC = g++ -m64

# Compiler flags
FLAGS = -O3 -Wall -g -fopenmp

# EXE file
EXE = $(PRJ_NAME) 

# Object files
OBJ = main.o power_grid.o sub_grid.o passive_elements.o current_source.o node.o interconnect_tree.o layer.o c4_pad.o profiler.o sparse_matrix.o

# Header files
H = ../include/power_grid.h ../include/profiler.h ../include/sparse_matrix.h

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) -o $(EXE) $(LIBPATHS) $(LIBS) $(LDOPT)

#$(PRJ_NAME).o: ../src/$(PRJ_NAME).cxx $(H)
#	$(CC) -c $(FLAGS) $(INCPATHS) ../src/$(PRJ_NAME).cxx

main.o: ../main.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../main.cxx

power_grid.o: ../src/power_grid.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/power_grid.cxx

sub_grid.o: ../src/sub_grid.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/sub_grid.cxx

passive_elements.o: ../src/passive_elements.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/passive_elements.cxx
	
current_source.o: ../src/current_source.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/current_source.cxx

node.o: ../src/node.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/node.cxx

interconnect_tree.o: ../src/interconnect_tree.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/interconnect_tree.cxx

layer.o: ../src/layer.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/layer.cxx

c4_pad.o: ../src/c4_pad.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/c4_pad.cxx

profiler.o: ../src/profiler.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) -I/usr/include/boost ../src/profiler.cxx

sparse_matrix.o: ../src/sparse_matrix.cxx $(H)
	$(CC) -c $(FLAGS) $(INCPATHS) ../src/sparse_matrix.cxx

#common.o: ../src/common.cxx $(H)
#	$(CC) -c $(FLAGS) $(INCPATHS) ../src/common.cxx

clean:
	rm -f *.o $(PRJ_NAME) 
