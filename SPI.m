clear
clc
close all
%%
set_calibration = 1;
set_1Dstitch = 1;
set_crop_rescale = 1;
set_Raster = 1;
set_upscale = 1;
set_WB = 1;

Default_path = "Y:\jia-lab\Kyungduck\Projects\TDIOSM\Presentations\SupplementarySoftware";

mkdir('Sample_data_1Dstitch')
mkdir('Instant_SPI_image')
mkdir('Raster_SPI_image')
mkdir('Upscaled_SPI_image')
mkdir('Final_SPI_image')

%% This software requires MIJ/MIJI plugin.
% mij.jar can be found at: https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab
% https://imagej.net/plugins/miji
% mij/ij files can be found in this folder
% This code use Grid/Collection stitching plugin


javaaddpath 'C:\Program Files\MATLAB\R2020b\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2020b\java\ij.jar'
% addpath('C:\Users\kyoon47\Downloads\fiji-win64\Fiji.app\scripts');

%% 0. Calibration parameters. Surface imperfections of MLA imposes FPN
cd(Default_path)
if set_calibration == 1
clc
disp('Setting up calibration...')
FilePath_calibration = '.\\Calibration';

offset = imread(strcat(FilePath_calibration,'\\AVG_4x1_I0.tif'));
I1 = imread(strcat(FilePath_calibration,'\\AVG_4x1_I1.tif'));

range=1301:2900; %% Check the X range of 2nd MLA from I1

FPN_profile = mean(I1-offset);

FPN_profile_range = FPN_profile(range);
Correction_weight = max(FPN_profile_range)./FPN_profile_range;
end
%% 1. Stitch buffers 
cd(Default_path)
if set_1Dstitch == 1
clc
disp('Stitching image buffers in 1D...')
data_struct = dir('.\Sample_data');
data_struct2 = data_struct(3:end);
folder_count = max(size(data_struct2))

% Call ImageJ
% MIJ.start
Miji();
% folder_directory = pwd;
folder_directory = 'Y:\\jia-lab\\Kyungduck\\Projects\\TDIOSM\\Presentations\\SupplementarySoftware\\Sample_data\\';
savepath_1D = 'Y:\\jia-lab\\Kyungduck\\Projects\\TDIOSM\\Presentations\\SupplementarySoftware\\Sample_data_1Dstitch\\';

gridX = 1;
gridY = folder_count;
overlap = 0;

for FOV = 1:gridY % Stitch along 1D scan
    
    x={dir(sprintf(strcat(folder_directory,data_struct2(FOV).name,'\\*.tiff'),FOV)).name};
    totalimgs=max(size(x));


    MIJ.run("Grid/Collection stitching", sprintf("type=[Grid: column-by-column] order=[Down & Left] grid_size_x=%d grid_size_y=%d tile_overlap=%d first_file_index_i=0 directory=Y:/jia-lab/Kyungduck/Projects/TDIOSM/Presentations/SupplementarySoftware/Sample_data/smear_%02d_4x1 file_names=TDI_1024.{ii}.tiff output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 display_fusion computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display] use",gridX,totalimgs-1,overlap,FOV-1));
    MIJ.run("Tiff...", sprintf(strcat("path=",savepath_1D,"Fused_1D_%02d.tiff"),FOV))
    MIJ.closeAllWindows
    clc
end
clc
disp('Stitching image buffers complete!')
end
%% 2. Crop, rescale 1D scan images
cd(Default_path)
if set_crop_rescale == 1
clc
disp('FPN correction...')
fname = ".\\Sample_data_1Dstitch\\Fused_1D_%02d.tiff";

SavePath = '.\\Instant_SPI_image';
sname = strcat(SavePath,'\\InstantSPI_%02d.tif');

for FOV = 1:gridY
    I = imread(sprintf(fname,FOV));
    I_cropped = I(:,range);
    I_offset = (I_cropped)-offset(range);
    I_flatfield = double(I_offset) .* Correction_weight/2;
    imwrite(uint16(I_flatfield),sprintf(sname,FOV),'tif');
end
clc
disp('FPN correction complete!')
end
%% 3. (Optional) Raster stitching
cd(Default_path)
if set_Raster == 1
clc
disp('Raster stitching computing 10% overlap')
savepath_Raster = 'Y:\\jia-lab\\Kyungduck\\Projects\\TDIOSM\\Presentations\\SupplementarySoftware\\Raster_SPI_image\\';
RasterX = folder_count;
RasterY = 1;
Raster_overlap = 10; % percentages
MIJ.run("Grid/Collection stitching", sprintf("type=[Grid: row-by-row] order=[Left & Down] grid_size_x=%d grid_size_y=%d tile_overlap=%d first_file_index_i=1 directory=Y:/jia-lab/Kyungduck/Projects/TDIOSM/Presentations/SupplementarySoftware/Instant_SPI_image file_names=InstantSPI_{ii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=3.50 compute_overlap absolute_displacement_threshold=3.50 display_fusion computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display] use",RasterX,RasterY,Raster_overlap,FOV-1));
MIJ.run("Tiff...", strcat("path=",savepath_Raster,"Raster_SPI.tif"))
MIJ.closeAllWindows
clc
disp('Stitching complete!')
end

%% 4. Upscaling
cd(Default_path)
% due to the vast image size exceeding its optimal operable range, 
% upscaling & deconvolution is demonstrated with 1D scan
if set_upscale == 1
    clc
    disp('Upscaling images...')
    savepath_Upscaled = 'Y:\\jia-lab\\Kyungduck\\Projects\\TDIOSM\\Presentations\\SupplementarySoftware\\Upscaled_SPI_image\\';
    imp = ij.ImagePlus("Y:/jia-lab/Kyungduck/Projects/TDIOSM/Presentations/SupplementarySoftware/Instant_SPI_image/InstantSPI_01.tif");
    imp.show()
    im = MIJ.getCurrentImage;
    [size_x, size_y] = size(im);
    MIJ.run("Scale...", sprintf("x=2 y=2 width=%d height=%d interpolation=Bilinear average create",2*size_x, 2*size_y));
    MIJ.run("Tiff...", strcat("path=",savepath_Upscaled,"SPI_upscaled.tif"))
    MIJ.closeAllWindows
    clc
    disp('Upscale complete!')
end

%% 5. Deconvolution
% Deconvolution 
% For WB Deconvolution, please use the method published by Guo et al. in
% Nature Biotechnology
% https://www.nature.com/articles/s41587-020-0560-x
% https://github.com/eguomin/regDeconProject/tree/master/WBDeconvolution

if set_WB ==1

data_path = 'Y:\jia-lab\Kyungduck\Projects\TDIOSM\Presentations\SupplementarySoftware\Upscaled_SPI_image\';
file_name = 'SPI_upscaled.tif';
psf_path = 'Y:\jia-lab\Kyungduck\PSF-37.tif';


full_name = strsplit(file_name,'.tif');
save_name = 'Final_SPI_image';
% set parameters
deconMethod = 2;
proMode = 0;
itNum = 1; % iteration number
t1 = 0;
t2 = 0;
gpuFlag = 0;
if(proMode==1)
    gpuFlag = 1;% 0: CPU; 1: GPU  
    gpuDevNum = 1; % specify the GPU device if there are multiple GPUs
end
path_output = ['Y:\jia-lab\Kyungduck\Projects\TDIOSM\Presentations\SupplementarySoftware\Final_SPI_image\'];
mkdir(path_output);
%%%%%%%%%%%%%%%%%%%%%%%% read in images %%%%%%%%%%%%%%%%%%%%%
% stackIn = single(ReadTifStack([data_path file_name]));
stackIn = single(ReadTifStack([data_path file_name]));
[Sx, Sy, Sz] = size(stackIn);

disp('Preprocessing forward and back projectors ...');
% % forward projector: PSF
PSFIn = single (ReadTifStack(psf_path));
PSF1 = PSFIn/sum(PSFIn(:));
% % % back projector: PSF_bp
% parameters: light sheet microscopy as an example
switch(deconMethod)
    case 1
        bp_type = 'traditional'; 
    case 2
        bp_type = 'wiener-butterworth';
    otherwise
        error('Processing terminated, please set deconconvolution method as 1 or 2')
end
alpha = 0.05;
beta = 1; 
n = 10;
resFlag = 1;
iRes = [2.44,2.44,10];
verboseFlag = 0;
[PSF2, ~] = BackProjector(PSF1, bp_type, alpha, beta, n, resFlag, iRes, verboseFlag);
PSF2 = PSF2/sum(PSF2(:));
WriteTifStack(PSF1, [path_output, 'PSF_fp.tif'], 32);
WriteTifStack(PSF2, [path_output, 'PSF_bp.tif'], 32);

% set initialization of the deconvolution
flagConstInitial = 0; % 1: constant mean; 0: input image

% % % deconvolution
PSF_fp = align_size(PSF1, Sx,Sy,Sz);
PSF_bp = align_size(PSF2, Sx,Sy,Sz);
if(gpuFlag)
    g = gpuDevice(gpuDevNum); 
    reset(g); wait(g); 
    disp(['GPU free memory: ', num2str(g.FreeMemory/1024/1024), ' MB']);
end
if(gpuFlag)
    OTF_fp = fftn(ifftshift(gpuArray(single(PSF_fp))));
    OTF_bp = fftn(ifftshift(gpuArray(single(PSF_bp))));
else
    OTF_fp = fftn(ifftshift(PSF_fp));
    OTF_bp = fftn(ifftshift(PSF_bp));
end
disp('Start deconvolution...');
smallValue = 0.001;
for imgNum = 1
    tic
    disp(['...Processing image #: ' num2str(imgNum)]);
    fileIn = [data_path file_name]; %
    stackIn = single(ReadTifStack(fileIn));
    if(gpuFlag)
        stack = gpuArray(single(stackIn));
    else
        stack = stackIn;
    end
    stack = max(stack,smallValue);
    if(flagConstInitial==1)
        stackEstimate = ones(Sx, Sy, Sz)*mean(stack(:)); % constant initialization
    else
        stackEstimate = stack; % Measured image as initialization
    end
    
    for i = 1:itNum
        stackEstimate = stackEstimate.*ConvFFT3_S(stack./...
        ConvFFT3_S(stackEstimate, OTF_fp),OTF_bp);
        stackEstimate = max(stackEstimate,smallValue);
    end
    if(gpuFlag)
        output = gather(stackEstimate);
    else
        output = stackEstimate;
    end
    toc
end
WriteTifStack(output, [path_output, sprintf(save_name), '.tif'], 16);

if(gpuFlag) % reset GPU
    reset(g); 
end
disp('Deconvolution completed !!!');

end

%%%%%%%%%%%%%%%%%%%%%%%%
% % % Function
function img2 = align_size(img1,Sx2,Sy2,Sz2,padValue)
if(nargin == 4)
    padValue = 0;
end

[Sx1,Sy1,Sz1] = size(img1);
Sx = max(Sx1,Sx2);
Sy = max(Sy1,Sy2);
Sz = max(Sz1,Sz2);
imgTemp = ones(Sx,Sy,Sz)*padValue;

Sox = round((Sx-Sx1)/2)+1;
Soy = round((Sy-Sy1)/2)+1;
Soz = round((Sz-Sz1)/2)+1;
imgTemp(Sox:Sox+Sx1-1,Soy:Soy+Sy1-1,Soz:Soz+Sz1-1) = img1;


Sox = round((Sx-Sx2)/2)+1;
Soy = round((Sy-Sy2)/2)+1;
Soz = round((Sz-Sz2)/2)+1;
img2 = imgTemp(Sox:Sox+Sx2-1,Soy:Soy+Sy2-1,Soz:Soz+Sz2-1);
end