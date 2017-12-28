TARGET = main
CC = g++-4.9 -m64
CFLAGS = -O3 -g -fopenmp -Wall 

# Mosek path
MOSEKPATH   = /autofs/fs1.ece/fs1.eecg.najm/b/b1/power_grid_code/mosek
MOSEKHEADER = $(MOSEKPATH)/8/tools/platform/linux64x86/h
MOSEKLIB    = $(MOSEKPATH)/8/tools/platform/linux64x86/bin
LDOPT = -Wl,-rpath-link,$(MOSEKPATH)/8/tools/platform/linux64x86/bin -Wl,-rpath,'/autofs/fs1.ece/fs1.eecg.najm/b/b1/power_grid_code/mosek/8/tools/platform/linux64x86/bin'
INCLUDE = inc/

# Destination of suitsparse libarires
LIBDEST = /autofs/fs1.ece/fs1.eecg.najm/b/b1/power_grid_code/libraries/lib

# Destination of suitsparse header files)
INCDEST = /autofs/fs1.ece/fs1.eecg.najm/b/b1/power_grid_code/libraries/include

# Link libraries
LIBS =  -lc -lm -lrt \
		$(MOSEKLIB)/libmosek64.so \
		$(LIBDEST)/libumfpack.a \
		$(LIBDEST)/libcholmod.a \
		$(LIBDEST)/libcolamd.a \
		$(LIBDEST)/libamd.a \
		$(LIBDEST)/libcsparse.a \
		$(LIBDEST)/libcxsparse.a \
		$(LIBDEST)/liblapack.a \
	    $(LIBDEST)/libsuitesparseconfig.a \
		$(LIBDEST)/libopenblas.a \
		$(LIBDEST)/libgfortran.a   


.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst src/%.cxx, obj/%.o, $(wildcard src/*.cxx))

obj/%.o: src/%.cxx 
	$(CC) $(CFLAGS) -c -I$(INCDEST) -I$(INCLUDE) -I$(MOSEKHEADER) $< -o $@ 

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $@ -I$(INCDEST) -I$(INCLUDE) -I$(MOSEKHEADER) $(LIBS) $(LDOPT)

clean:
	-rm -f obj/*.o $(TARGET)
