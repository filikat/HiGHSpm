## Notes on installing spral - for SSID

 - Needs hwloc from Git or WWW, and that seems to need pkgconfig

 - Needs gklib from Git

 - Needs METIS from Git

Get spral from Git

ensure that autoreconf is available (sudo aopt install autoreconf)

May need

 - apt-get install libblas-dev
 - apt-get install liblapack-dev

./autogen

spral_complex.h needs to be added to /usr/local/include/ by hand

When running code that calls SSID, must have

 - export OMP_CANCELLATION=true
 - export OMP_PROC_BIND=true
