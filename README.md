Module for reading out and plotting data from MESA binary simulations, intended for use in tandem with MESA Binary from iadd.astro.illinois.edu. 
Most of the functions produce animations (at the time of writing, only the kippenhahn plot is a static plot). 
Users can select the desired FPS of the movie (must be an int) and for specific plots (like abundances), threshold values must be specified (also must be an int).

For Personal Use:

Given the source code, the methods used here could easily be modified for single star evolution. Simply remove the "logs2" from the functions and second line calling the plotting functions. 

To streamline the process of aligning timesteps, this code assumes that history_interval and profile_interval are equal for star 1, star 2, and the binary data. 
On iADD's MESA Binary, this is done by default. If you want to use this code for simulations run on your own machine, make sure history_interval and profile_interval are equal for both stars and the binary data. 

For the kippenhahn plot, you must edit profile_columns.list to include gradr (radiative gradient) and grada (adiabatic gradient). This is done by simply removing the ! before their entries.
For the abundance plot to work, you must edit profiles_columns.list to include the elements/isotopes you want. Again, this is done by simply removing the ! before their entries.

For questions, comments, and concerns about this code, please contact njhurst2@illinois.edu
