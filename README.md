# ZymographyCode
Methodology for image segmentation proposed in the article "Quantification of the presence of enzymes in gelatin zymography using the Gini Index".

You must have installed the free software R. It can be downloaded at the following link:
https://www.r-project.org/

A) Once the software is installed, open the PACKAGES.R file to install the necessary packages for the program.
   This step could take a few minutes.

B) Run the functions of the following files: (IMPORTANT: Ctrl + Alt + R – Runs the entire script)

	1- MIGB.R - GI algorithm function to adjust "s" Gaussians in 2-dimensional data set
  
	2- Background.R - Function for background analysis
  
	2- UmB.R  - Functions to obtain the threshold value "Um" and the "B" value
  
	3- BasicFunctions - Functions used to perform the steps of the proposed methodology

C) Open the MainCode.R file and follow the instructions to use the proposed methodology.
