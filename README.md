# gavrikov-detonation-cell-size-model
This is an implementation of Gavrikov's et al. detonation cell size model in Cantera. Source article:
Gavrikov A. I., Efimenko A. A., and Dorofeev S. B., A Model for Detonation Cell Size Prediction from Chemical Kinetics, Russian Research Center, Kurchatov Institute, Russia, 1999

Exemplary usage is shown in "example.py" scipt, which calculates detonation cell size for a stoichiometric hydrogen-air mixture.

The three main steps of calculations are as follows:
* calculate the reaction lengths defined at 0.75 Ma using the ZND model, as well as the von Neumann temperature
* calculate the effective activation energies as defined in Gavrikov's work along with the post shock temperature
* calculate the detonation cell size using Gavrikov's formula

In order to run this script, you need to have the newest Cantera and SDToolbox already installed. More information under the links below.

Cantera:
http://www.cantera.org/docs/sphinx/html/index.html

SDToolbox:
http://shepherd.caltech.edu/EDL/public/cantera/html/SD_Toolbox/

TODO:
* better ode integration
