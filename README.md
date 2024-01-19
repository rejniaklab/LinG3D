# Project Title

## LinG3D: Visualizing the Spatio-Temporal Dynamics of Clonal Evolution

## Description

This code generates the 3D lineage trees of (1) all clones; (2) individual clones; (3) all clones , but with only those cells that survived to the end of simulation; and (4) individual clones containing only those cells that survived to the end of simulation – see the quick guide to LinG3D routines. The code is generated in three programming languages: Matlab, R, and Python.

## Dependencies

MATLAB software (the R2020b version on a Mac computer was used for all testing)

R language (the R version 4.1.1 on a Mac computer (macOS Big Sur) was used for all testing)

Python (version >3.8 on a Mac computer was used for testing)

## Installing

Matlab: this code does not need installation

R: the required libraries are readr, rapportools, rgl, and devtools. 

        install.packages("readr")   

        install.packages("rapportools")
        
        install.packages("rgl")

        install.packages("devtools")
        
Python: It requires the following libraries:

         Numpy
        
         Scipy
         
         Matplotlib
         
## Executing the program**

Matlab: 
Download all MATLAB files (add all files and folders to your MATLAB path)
Run one of the routines listed below. 

R:

OPTION 1 (Install the LinG3D package)
   
        library(devtools)  # load the devtools library

        install_github("rejniaklab/r_LinG3D")  # install the package from GitHub

        library(LinG3D)  # load the LinG3D library

        linG3DClone(arguments)  # run the routine (example)

OPTION 2 (Use the routines)

        1.	Download all R codes and data folders to your computer.
        
        2.	Source each function. For example, to source linG3DClone.R saved in your_folder on the desktop type on the console,
                source("~/Desktop/your_folder/linG3DClone.R").
                
        3.	Run the function. For example, to run linG3DClone.R type on the console,
                linG3DClone(arguments).
                          
Python:
       You can run the .py scripts from the terminal
       
       Terminal> git clone https://github.com/rejniaklab/LinG3D.git

       Terminal> cd LinG3D

       Terminal> python3 linG3DAliveAll.py
       
       Terminal> python3 linG3DAliveClone.py
       
       Terminal> python3 linG3DAll.py
       
       Terminal> python3 linG3DClone.py
     
       The codes are also available as Jupyter notebooks


## File manifest

The following set of functions are available, and explained in detail in the quick guide to LinG3D routines.
i)	linG3DAliveAll
ii)	linG3DAliveClone
iii)	linG3DAll
iv)	linG3DClone

## Example
    
```python
linG3DClone(pathData='exampleB05',cloneNum=5,IsGradient = 1,xmin=-100,xmax=100,ymin=-100,ymax=100,tmin=0,tmax=100000,fileStep = 2000,toPrint=1)
```

<div style="margin:2%";>  
    <img src="https://github.com/okayode/pyLinG3D/blob/okayode/exampleB05/fig_clones/tree_clone_5.jpg?raw=true"; alt="tree_clone_5"; width=30%;/>
</div>

## Authors

Anjun Hu,
Awino Maureiq E. Ojwang’,
Kayode Olumoyin,
Katarzyna Rejniak

## Version
0.1
Initial release

## License
This project is licensed under the GNU General Public License v3.0.
