# SPI Image Reconstruction
This project focuses on the reconstruction pipeline of images obtained from the Super-resolution Panoramic Integration (SPI) Microscopy. A main SPI.m source code contains all necessary data processing included in the project.
The code follows five main steps: 
1. 1D stitching of swept fields of view (FOV)
2. Image crop & FPN removal
3. (Optional) Large FOV formation for raster-scanned images
4. Upscaling
5. WB Deconvolution


# Prerequisites:
This pipeline require MIJ, a Java package for running ImageJ and Fiji within Matlab. While a relevant files are uploaded within the package, they could be found at:
https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab
https://imagej.net/plugins/miji
Grid/Collection stitching plugin by Preibisch et al., Bioinformatics (2009) is used in this source code.
Prior to running the SPI.m source code, please add mij.jar and ij.jar to the java folder of your MATLAB path. 
e.g. 
javaaddpath 'C:\Program Files\MATLAB\R2020b\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2020b\java\ij.jar'

# Steps:
Step 1. The size of raw images collected by the TDI sensor is determined by the buffer height set by the frame grabber. As the SPI system images super-resolution samples in motion, images need to be stiched to represent a true specimen captured by sample scanning. Example raw images are located in Sample_data folder.

Step 2. Since the width of the TDI sensor (4k) vastly surpasses the microscope FOV, an additional step is required to crop out unnecessary FOV captured by the sensor. Furthermore, due to a Gaussian collimated beam illuminating the MLA, optical post processing results in weaker signals near the periphery and surface imperfections further impose fixed pattern noise (FPN) on the instant image collected by SPI. Sample calibration images are provided in the Calibration folder, time-averaged images of dark signal and a slide of uniform fluorescent dye.

Step 3. In case of raster scanning, SPI can obtain an image of unconstrained FOV suitable for high-throughput imaging applications. We use Grid/Collection stitching plugin in ImageJ.

Step 4. Example images have 74-nm pixel resolution. While this is sufficient for sampling intermediate 1.4x improved resolution beyond the diffraction limit (Rayleigh criterion), it is not sufficient sampling a full twofold improved resolution. A simple image upscaling is implemented to be performed before deconvolution.

Step 5. WB deconvolution accelerates processing time by approximately 40 fold, which is suitable for high-throughput images obtained from SPI. We integrated WB Deconvolution functions into our source code. (https://github.com/eguomin/regDeconProject/tree/master/WBDeconvolution)

# Sample images:
Sample images can be found at:
https://gtvault-my.sharepoint.com/:f:/g/personal/kyoon47_gatech_edu/EtnweSat9TZBietUQ_cuK0kBo4RyfpOQikDBV-c41kIl8g?e=my6l4H

# References:
Please refer to relevant academic papers, documentation, or resources for a deeper undestanding of the SPI system and the reconstruction techniques used.
For questions or assistance, please contact Kyungduck Yoon at leo.k.yoon@gatech.edu
