# DaXi: High-Resolution, Large Imaging Volume, and Multi-View Single Objective Light-Sheet Microscope

![DaXi_blender_render](https://user-images.githubusercontent.com/1870994/223275335-062fcd03-7a1a-4cdb-8c76-c9e8312d6466.png)

The DaXi microscope is a single-objective light-sheet microscope design that combines the advantages of standard single-objective microscopes with the speed, coverage, resolution, and gentleness of light-sheet microscopes. Its main features include a wider field of view and high-resolution imaging via a custom remote focusing objective, fast volumetric imaging over larger volumes without compromising image quality or necessitating tiled acquisition, fuller image coverage for large samples via multi-view imaging, and higher throughput multi-well imaging via remote coverslip placement. The microscope achieves a resolution of 450 nm laterally and 2 μm axially over an imaging volume of 3,000 × 800 × 300 μm, making it suitable for imaging large living samples such as developing embryos. The DaXi microscope also improves sample mounting ergonomics and can be converted from upright to inverted using remote focusing. The microscope has been successfully used to image various systems, including Drosophila egg chamber development, zebrafish whole-brain activity, and zebrafish embryonic development, up to nine embryos at a time.

DaXi has been published in Nature Methods in 2022:
[Yang, B., Lange, M., Millett-Sikking, A. et al. DaXi—high-resolution, large imaging volume and multi-view single-objective light-sheet microscopy. Nat Methods 19, 461–469 (2022).](https://doi.org/10.1038/s41592-022-01417-2)

Biorxiv preprint: https://doi.org/10.1101/2020.09.22.309229
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

Bin Yang1*, Merlin Lange, Alfred Millett-Sikking, Ahmet Can Solak, Shruthi Vijaykumar, Wanpeng Wang, 
Hirofumi Kobayashi, Matthew N Mccarroll, Lachlan W Whitehead, Reto P Fiolka, Thomas B. Kornberg, 
Andrew G York, Loic A. Royer*

Link to the Biorxiv preprint:  https://doi.org/10.1101/2020.09.22.309229


## License
This project is covered under the BSD 3-clause license.

Copyright © 2021. Chan Zuckerberg Biohub. All rights reserved.

![image](https://user-images.githubusercontent.com/1870994/115829677-1c617a80-a3c4-11eb-9011-200e3e0ab435.png)
