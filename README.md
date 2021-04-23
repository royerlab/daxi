![image](https://user-images.githubusercontent.com/1870994/115829677-1c617a80-a3c4-11eb-9011-200e3e0ab435.png)

# DaXi: High-Resolution, Large Imaging Volume, and Multi-View Single Objective Light-Sheet Microscopy

Recent developments in Oblique Plane Microscopy (OPM) have shown that it can achieve high spatio-temporal resolution. Here we describe a single objective light-sheet microscope based on oblique plane illumination that achieves: (i) large field of view and high-resolution imaging via a custom remote focusing objective; (ii) fast volumetric imaging by means of light-sheet stabilised stage scanning – a novel scanning modality that extends the imaging volume without compromising imaging speed nor quality; (iii) multi-view imaging by alternating the orientation of light-sheet illumination and detection to improve the image quality on large samples; (iv) simpler design and ergonomics by remote placement of coverslips to allow inverted imaging, enabling imaging across scales in a high-throughput format. Overall, we achieved a resolution of 450 nm laterally and 2 μm axially and a field of view of 3000 μm × 800 μm × 300 μm. We demonstrate the speed, field of view, resolution and versatility of our novel instrument by imaging various systems, including zebrafish whole brain activity, Drosophila egg chamber development, and zebrafish development – up to nine embryos simultaneously.

Link to the Biorxiv preprint: https://doi.org/10.1101/2020.09.22.309229
Suplementary information for the paper can be found [here](https://www.biorxiv.org/content/10.1101/2020.09.22.309229v2.supplementary-material).  

This repository contains all CAD drawings, parts list, control software, and simulation for the microscope.
It is our goal to enable the community to build their own replicates of this system. Don't hesitate to post 
issues to this repo to request addtional information or report issues.  

Correspondence: loic.royer@czbiohub.org, bin.yang@czbiohub.org

## Contents

* [CAD models](https://github.com/royerlab/daxi/tree/master/cad_models): custom designed parts used to built the microscope. The STL files can be used directly 
for 3D printing.

* [Part list](https://github.com/royerlab/daxi/tree/master/part_list): a list of all the necessary parts used to built the microscope. 
A rough estimate of the cost of the microscope is about $130k.

* [Control code](https://github.com/royerlab/daxi/tree/master/control_code): to control the microscope. A Micro-Manager script to setup the data acquisition and 
a python module to programme a NI DAQ system to synchronize other devices during the acquisition. 
Micro-Manager is used to control the camera and the microscope stage. 
The NI DAQ system uses digital signal to switch the lasers on/off, and analog signals to set the angle of
the galvo mirrors and the position of the piezo actuators. 

* [Simulation code](https://github.com/royerlab/daxi/tree/master/simulation_code): used to generate the results in Supplementary Figures 7 and 9. 

## Image Processing

The image processing code required to go from DaXi's raw data to a fully processed timelapse or even video is found [here](https://github.com/royerlab/dexp).

## Software Requirements
Code is only tested on Windows 10 operating system, many hardware drivers only exist for windows...

Micro-Manager Version is 2.0-gamma.

Python version 3.7.

## Citation:

If you use the hardwar designs or software code, please cite us:

[High-Resolution, Large Imaging Volume, and Multi-View Single Objective Light-Sheet Microscopy](https://doi.org/10.1101/2020.09.22.309229)

Bin Yang1*, Merlin Lange1, Alfred Millett-Sikking2, Ahmet Can Solak1, Shruthi Vijaykumar1, Wanpeng Wang3, 
Hirofumi Kobayashi1, Matthew N Mccarroll4, Lachlan W Whitehead5;6, Reto P Fiolka7;8, Thomas B. Kornberg3, 
Andrew G York2, Loic A. Royer1*

1. Chan Zuckerberg Biohub, San Francisco, USA
2. Calico Life Sciences LLC, South San Francisco, USA
3. Cardiovascular Research Institute, University of California, San Francisco, CA 94143
4. Institute for Neurodegenerative Diseases, University of California, San Francisco, CA, 94143, USA.
5. The Walter and Eliza Hall Institute of Medical Research, Parkville, VIC, Australia
6. Department of Medical Biology, The University of Melbourne, Parkville, VIC, Australia
7. Department of Cell Biology, University of Texas Southwestern Medical Center, Dallas, United States
8. Lyda Hill Department of Bioinformatics, University of Texas Southwestern Medical Center, Dallas, United States

Link to the Biorxiv preprint:  https://doi.org/10.1101/2020.09.22.309229



## License
This project is covered under the BSD 3-clause license.

Copyright © 2021. Chan Zuckerberg Biohub. All rights reserved.
