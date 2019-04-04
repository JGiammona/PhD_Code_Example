# Jamming Simulator makefile

# Load the common configuration file
include config.mk

#SRC=simdefs.cc
#OBJ=$(SRC:.cc=.o)

# List of executables
EXECUTABLES= C_elegans_simulator

# Makefile rules
all: $(EXECUTABLES)

C_elegans_simulator: C_elegans_simulator.o printfuncs.o simdefs.o update.o
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) $^ -o $@ -lvoro++

C_elegans_simulator.o: C_elegans_simulator.cc simdefs.hh printfuncs.hh update.hh
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -c $< -o $@


printfuncs.o: printfuncs.cc printfuncs.hh simdefs.hh
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -c $< -o $@

simdefs.o: simdefs.cc simdefs.hh
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -c $< -o $@

update.o: update.cc simdefs.hh update.hh printfuncs.hh
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -c $< -o $@
clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
