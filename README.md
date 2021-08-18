# compositeHeatConduction
How to run: 
output files: *.temp.
plotting : run python3 animation.py which will generate a movie file animation.mp4.

The heat conduction code is written in object oriented way. Each shell is an object.
A shell may by rectangular or axisymmetric. By default it will be axisymmetric (a cylindrical shell). The can be changed by the member function void cylindrical(bool)
Properties of the shell are Length, Width, inner and outer radius (only for cyl. shell), K, rho, Cp etc.. and can be changed by set member functions.
Option for variable thermal conductivity w.r.t temperature is included.
Thermal conductivity along both the direction are independent. tk_x and tk_y are defined.
Intershell bc are : Contact resistance, Radiation, Radiation + Convection.

Shell system is defined like matrix elements. Code automatically updated length column wise, taking the value of the first shell in the column. likewise, do for width 
along rows. 
Simulation time is a shell property. For shell system, it copies from the first row, first column shell.
