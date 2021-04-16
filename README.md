# daxi supplementary information

Supplementary information related to the paper reporting daxi, a single-objective light
sheet microscopy method. 

## About the paper

High-Resolution, Large Imaging Volume, and Multi-View Single Objective Light-Sheet Microscopy

Bin Yang1*, Merlin Lange1, Alfred Millett-Sikking2, Ahmet Can Solak1, Shruthi Vijaykumar1, Wanpeng Wang3, Hirofumi Kobayashi1, Matthew N Mccarroll4, Lachlan W Whitehead5;6, Reto P Fiolka7;8, Thomas B. Kornberg3, Andrew G York2, Loic A. Royer1*

1 Chan Zuckerberg Biohub, San Francisco, USA

2 Calico Life Sciences LLC, South San Francisco, USA

3 Cardiovascular Research Institute, University of California, San Francisco, CA 94143

4 Institute for Neurodegenerative Diseases, University of California, San Francisco, CA, 94143, USA.

5 The Walter and Eliza Hall Institute of Medical Research, Parkville, VIC, Australia

6 Department of Medical Biology, The University of Melbourne, Parkville, VIC, Australia

7 Department of Cell Biology, University of Texas Southwestern Medical Center, Dallas, United States

8 Lyda Hill Department of Bioinformatics, University of Texas Southwestern Medical Center, Dallas, United States

Correspondence: loic.royer@czbiohub.org, bin.yang@czbiohub.org

All supplementary information related to the paper can be found here.

Link to the Biorxiv preprint: https://www.biorxiv.org/content/10.1101/2020.09.22.309229v1


## Repo Contents

* Cad models: custom designed parts used to built the microscope. The STL files can be used directly 
for 3D printing.

* Part list: a list of all the necessary parts used to built the microscope

* Control code: to control the microscope. A Micro-Manager script to setup the data acquisition and 
a python module to programme a NI DAQ system to synchronize other devices during the acquisition. 
Micro-Manager is used to control the camera and the microscope stage. 
The NI DAQ system uses digital signal to switch the lasers on/off, and analog signals to set the angle of
the galvo mirrors and the position of the piezo actuators. 

* Simulation code: used to generate the results in Supplementary Figures 7 and 9. 


## Software Requirements
Code is only tested on Windows 10 operating system. 

Micro-Manager Version is 2.0-gamma.

Ptyhon version 3.7.