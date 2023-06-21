
Before compiling the code, open script.sh and edit the variables:

- YOUR_CPP_COMPILER with the c++ compiler available
- YOUR_PATH_TO_INCLUDE with the path to reach include/highs/
- YOUR_PATH_TO_LIBRARY with the path to reach lib/libhighs

For Mac OS, it may be required to link other libraries (in AUX_LIBRARY_FOR_MAC). 
Otherwise, remove the line  "-L$AUX_LIBRARY_FOR_MAC -lz".

Add any compiler flag to FLAGS.

Then, type "bash script.sh" in a terminal. This creates the executable "ipm".

To solve an LP stored in ###.mps, type "./ipm ###.mps".

Example: 
    bash script.sh
    ./ipm afiro.mps
    ./ipm adlittle.mps

