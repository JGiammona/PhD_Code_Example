# PhD_Code_Example

**Example_code** shows an example of the code I wrote during my PhD to examine the robustness of the early _C._ _elegans_ embryo to noise in division timing, division angles and offsets on when different cells divide

_C_elegans_simulator.cc_ is the main simulator code.

_update.cc_ contains the definitions of many of the simulator functions.

_simdefs.cc_ contains some miscellaneous function definitions.

_printfuncs.cc_ contains printing function definitions.

Run "make" in this directory to compile code. 

Run "./C_elegans_simulator" to run.

Running "render timeseries.nb" with the Import[] functions adjusted to your machine's directory path will show a rendering of your simulation in **C_elegans_test**
