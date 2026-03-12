# iADD MESA Binary (WIP)

## Background

Package for reading out and plotting data from MESA binary simulations, designed in tandem with MESA Binary from [iadd.astro.illinois.edu](https://iadd.astro.illinois.edu/). Users are intended to run MESA binary simulations and then create plots as desired with this Package. Please see [this page](https://iadd.astro.illinois.edu/mesa-binary-background) for more information and the MESA binary simulator.

Most of the functions produce animations (at the time of writing, only the kippenhahn plot is a static plot), giving users the ability to watch binary stellar evolution through multiple properties.
Users can select the desired FPS of the movie and for specific plots (like abundances), threshold values must be specified. Movies are designed to be formatted such that they are accessible to all; colorblind friendly, large font sizes, etc. However, if you wish to see anything added, please contact us at the information found at the bottom of this page. 

## Installation & Usage

To install this module, use pip:

```
pip install iadd-mesa-binary
```

## Additional Dependencies

Users will need to have the following installed for some or all of the plots. Clicking on the hyperlink takes you to their installation pages, most of which are simply a pip install command.

* [mesaPlot](https://github.com/rjfarmer/mesaplot)
* [mesa_reader](https://github.com/wmwolf/py_mesa_reader)
* [TULIPS](https://astro-tulips.readthedocs.io/en/latest/installation.html) (needed for abundances plot)
* [Latex](https://www.tug.org/texlive/acquire-netinstall.html) (needed for text formatting, optional)

## Functions

Behind all of these functions are auxillary helper functions not listed here. They are however included in the download above, if you desire to look at them.

logs1 is the path to the logs folder for star 1. Must be a string.\
logs is the path to the logs folder for star 2. Must be a string.\
binary_history is the path to the binary_history.data file. Must be a string.\
fps is the desired frames per second of the generated movie. Generally, 8-16 fps makes for good movies. Must be an int.

```
plot_Mass_Transfer(logs1, logs2, binary_history, fps):
```

```
plot_Roche_Lobe(logs1, logs2, binary_history, fps):
```

```
plot_Hertzsprung_Russel(logs1, logs2, binary_history, fps, observers):
```
Observers is a boolean; set to True, it will slightly modify the HR diagram to include spectral classes and a gradiented colorbar at the bottom.

```
plot_Abundances(logs1, logs2, binary_history, fps, threshold):
```
If an element is below the given treshold value, it is removed from the plot. 0.001 is typically very reasonable. Must be given as an int. 

```
plot_Kippenhahn(logs1, logs2):
```
Unlike the rest of the above functions, plot_Kippenhahn produces a picture, not a movie. 

## Personal Use

Given the source code, the methods used here could easily be modified for single star evolution. Simply remove the "logs2" from the functions and second line calling the plotting functions. 

To streamline the process of aligning timesteps, this code assumes that history_interval and profile_interval are equal for star 1, star 2, and the binary data. On iADD's MESA Binary, this is done by default. If you want to use this code for simulations run on your own machine, make sure history_interval and profile_interval are equal for both stars and the binary data. 

For the kippenhahn plot, you must edit profile_columns.list to include gradr (radiative gradient) and grada (adiabatic gradient). This is done by simply removing the ! before their entries.
For the abundance plot to work, you must edit profiles_columns.list to include the elements/isotopes you want. Again, this is done by simply removing the ! before their entries.

For questions, comments, and concerns about this code, please contact njhurst2@illinois.edu.



