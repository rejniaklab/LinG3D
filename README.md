**Project Title**
LinG3D: Visualizing the Spatio-Temporal Dynamics of Clonal Evolution

**Description**
This code generates the 3D lineage trees of (1) all clones; (2) individual clones; (3) all clones , but with only those cells that survived to the end of simulation; and (4) individual clones containing only those cells that survived to the end of simulation – see the quick guide to LinG3D routines. The code is generated in three programming languages: Matlab, R, and Python.

**Dependencies**

MATLAB software (the R2020b version on a Mac computer was used for all testing)

R language (the R version 4.1.1 on a Mac computer (macOS Big Sur) was used for all testing)

Python (version 3.10.6 on a Mac computer was used for testing)

**Installing**

Matlab: this code does not need installation

R: the required libraries are readr, rapportools, and rgl. 

        install.packages("readr")   

        install.packages("rapportools")
        
        install.packages("rgl")

        install.packages("devtools")
        
Python: It requires the following libraries:

         Numpy
        
         Scipy
         
         Matplotlib
         
**Executing from R package**

library(devtools)

install_github("rejniaklab/r_LinG3D")

library(LinG3D)

**Executing the program**

Matlab: 
Download all MATLAB files (add all files and folders to your MATLAB path)
Run one of the routines listed below. 

R:
1.	Download all R codes data folders to your computer
2.	Source each function. For example, to source linG3DClone.r saved in your_folder on the desktop type on the console,
             source("~/Desktop/your_folder/linG3DClone.r")
3.	Run the function. For example, to run linG3DClone.r type on the console,
                          linG3DClone()
                          
Python:
       You can run the .py scripts from the terminal
       
       Terminal> python3 linG3DAliveAll.py
       
       Terminal> python3 linG3DAliveClone.py
       
       Terminal> python3 linG3DAll.py
       
       Terminal> python3 linG3DClone.py
     
       The codes are also available as Jupyter notebooks


**File manifest**
The following set of functions are available, and explained in detail in the quick guide to LinG3D routines.
i)	linG3DAliveAll
ii)	linG3DAliveClone
iii)	linG3DAll
iv)	linG3DClone

**Authors**
Anjun Hu,
Maureiq Ojwang’,
Kayode Olumoyin,
Katarzyna Rejniak

**Version**
0.1
Initial release

**License**
This project is licensed under the GNU General Public License v3.0.
