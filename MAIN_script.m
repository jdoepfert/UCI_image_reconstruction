%% Load data and set parameters

load data

Npix = 128;    % number of pixels
Np_full = 200; % number of projections for fully sampled data sets
Np_sub = 10;   % number of projections for subsampled data sets
Nimages = 15;  % number of images

addpath(genpath('NUFFT')); 
%% Reconstruct sub-sampled UCI data set (see Fig. 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Generate signal mask from S_off scan, for LS constraint                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
display('---------------------------------------------------------------------------')
display('Generating signal mask from Soff scan, this is needed for the LS constraint');
display('---------------------------------------------------------------------------')

% merge data from all slices into a single k-space
rawdata_S_off_sub_merged = ...
    reshape(rawdata_S_off_sub(:), Npix, Np_sub * Nimages);

% merge k-space trajectories from all slices  to a single trajectory
k_traj_UCI_merged = reshape(k_traj_UCI_sub, Npix, Np_sub * Nimages);

% plot trajectories
figure;
subplot(1,2,1)
plot(k_traj_UCI_sub(:,1),'b.','MarkerSize',4); axis image; axis off
title('subsampled trajectory of first slice')

subplot(1,2,2)
plot(k_traj_UCI_merged,'b.','MarkerSize',4); axis image; axis off
title('merged trajectory of all slices')

display('This is the subsampled k-space trajectory for a single slice. By merging the')
display('trajectory of all slices, a denser sampling is achieved.');
display('***Press button to continue***'); waitforbuttonpress;

% calculate density compensation
r = linspace(-0.5, 0.5, Npix)';
w = abs(r) * Npix / Np_sub;
w = repmat(w, [1, Np_sub*Nimages]);  

% calculate NUFFT of merged k-space
FT = NUFFT(k_traj_UCI_merged, w, 1, 0, [Npix, Npix], 2);
Fullimage = FT' * rawdata_S_off_sub_merged;

% create mask by thresholding
Fullimage_gray = mat2gray(abs(Fullimage));
mask_thresh = Fullimage_gray > 0.2; 

% only retain the three largest areas in mask with connectivity analysis
[L,num] = bwlabel(mask_thresh);   % determine connected segments
for i=1:num                       % calculate area of each segment
    area(i) = bwarea(L == i);
end
[~, idx] = sort(area,'descend');  % get indizes sorted regarding area size
I = idx(1:3);                     % get indizes of the three largest areas
mask = zeros(size(mask_thresh));
for i = 1:3                       % create mask from areas
    mask = mask+(L==I(i));        
end

% plot signal mask
subplot(1,2,1)
imagesc(Fullimage_gray); colormap gray; axis off; axis image
title('NUFFT image from merged k-space')
subplot(1,2,2)
imagesc(mask); colormap gray; axis off; axis image
title('Binary signal mask')

display('Therefore, the merged raw data can be well reconstructed by Fourier transform,')
display('leading to an acceptable signal mask.');
display('***Press button to continue***'); waitforbuttonpress;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Create imaging operator A and its adjoint At                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('---------------------------------------------------------------------------------')
display('Calculate standard Fourier transform reconstruction of subsampled UCI data set...');
display('---------------------------------------------------------------------------------')

% calculate density compensation
r = linspace(-0.5, 0.5, Npix)';
w = abs(r) * Npix / Np_sub;
w = repmat(w, [1, Np_sub]);

% create a NUFFT transform object for each image
for i = 1:Nimages
    NUFFT_objects{i} = ...
        NUFFT(reshape(k_traj_UCI_sub(:,i),Npix,Np_sub), w, 1, 0,[Npix,Npix],2);
end

% create imaging operators for entire data set
A = @(z)A_NUFFT(z, NUFFT_objects); 
At = @(z)At_NUFFT(z, NUFFT_objects);

% show simple NUFFT reconstruction
reco_NUFFT_S_on_sub = At(rawdata_S_on_sub);
reco_NUFFT_S_off_sub = At(rawdata_S_off_sub);

subplot(1,2,1)
ShowImages(abs(reco_NUFFT_S_on_sub));
title('FT reco of subsampled UCI S_o_n scan')

subplot(1,2,2)
ShowImages(abs(reco_NUFFT_S_off_sub));
title('FT reco of subsampled UCI S_o_f_f scan')

display('The Fourier transform reconstruction of the subsampled UCI data set exhibits');
display('strong undersampling artefacts.');
display('***Press button to continue***'); waitforbuttonpress;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Do SVT-LS reconstruction                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('----------------------------')
display('Do the SVT-LS reconstruction');
display('----------------------------')

% set parameters
k_max = 250;            % maximum number of iterations
delta = 0.8;            % regulates data fidelity
lambda = 0.76;          % regulates LS constraint
tau = 0.00002;          % regulates singular value shrinkage
tol = 0.09;             % relative tolerance level


% start reconstruction for S_on and S_off scan

display('Therefore, do SVT-LS reconstruction instead.');
display('Reconstructing Son...');
tic
    reco_S_on_sub = ...
        SVT_LS(rawdata_S_on_sub, A, At, tau, delta, lambda, k_max, tol, mask);
toc

display('Reconstructing Soff...');
tic
    reco_S_off_sub = ...
        SVT_LS(rawdata_S_off_sub, A, At, tau, delta, lambda, k_max, tol, mask);
toc

% plot reco
subplot(1,2,1)
ShowImages(abs(reco_S_on_sub));
title('SVT-LS reco of subsampled UCI S_o_n scan')

subplot(1,2,2)
ShowImages(abs(reco_S_off_sub));
title('SVT-LS reco of subsampled UCI S_o_f_f scan')

display('Most of the artefacts could be suppressed by the SVT-LS reconstruction.');
display('***Press button to continue***'); waitforbuttonpress;

%% Reconstruct full UCI data set for comparison (see Fig. 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Create imaging operator A and its adjoint At for full data set                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('----------------------------------------------------')
display('Compare subsampled UCI to fully sampled UCI data set');
display('----------------------------------------------------')
display('Reconstruct fully sampled UCI data set with Fourier transform...');

% calculate density compensation
r = linspace(-0.5, 0.5, Npix)';
w = abs(r) * Npix / Np_full;
w = repmat(w, [1, Np_full]);

% create a NUFFT transform object for each image
for i = 1:Nimages
    NUFFT_objects{i} = ...
        NUFFT(reshape(k_traj_UCI_full(:,i),Npix,Np_full), w, 1, 0,[Npix,Npix],2);
end

% create imaging operators for entire data set
A_full = @(z)A_NUFFT(z, NUFFT_objects); 
At_full = @(z)At_NUFFT(z, NUFFT_objects);

% show reconstruction
reco_NUFFT_S_on_full = At_full(rawdata_S_on_full);
reco_NUFFT_S_off_full = At_full(rawdata_S_off_full);

subplot(1,2,1)
ShowImages(abs(reco_NUFFT_S_on_full));
title('FT reco of full UCI S_o_n scan')

subplot(1,2,2)
ShowImages(abs(reco_NUFFT_S_off_full));
title('FT reco of full UCI S_o_f_f scan')

display('As expected, the fully sampled data sets show no subsampling artefacts.');
display('***Press button to continue***'); waitforbuttonpress;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Compare Son scan of full UCI to Son scan of subsampled UCI                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% selected images to compare
indizes = [1,3,4,8];
ref = abs(reco_NUFFT_S_on_full(:,:,indizes));
sub = abs(reco_S_on_sub(:,:,indizes));

% normalize images according to their maximum values
normalizer = max(ref(:)) / max(sub(:));

% show images
subplot(3,1,1)
ShowImages(ref, 1);
title('UCI S_o_n full')

subplot(3,1,2)
ShowImages(sub * normalizer, 1);
title('UCI S_o_n subsampled')
CLIM = get(gca, 'CLim');

subplot(3,1,3)
h = ShowImages(abs(ref - sub * normalizer),1);
title('difference')
set(gca,'CLim', CLIM); 

display('Comparing selected images from fully sampled and subsampled UCI demonstrates')
display('the good performance of the SVT-LS reconstruction');
display('***Press button to continue***'); waitforbuttonpress;

%% Compare CEST spectra (see Fig. 4)
% This piece of code requires "roipoly.m" from the Image Processing Toolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Reconstruct standard CEST data set                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('--------------------------')
display('Comparing the CEST spectra');
display('--------------------------')
display('Reconstruct standard CEST data set with Fourier transform...');

% calculate density compensation
r = linspace(-0.5, 0.5, Npix)';
w = abs(r) * Npix / Np_full;
w = repmat(w, [1, Np_full]);

% do FT
for i = 1:Nimages          
    FT = NUFFT(reshape(k_traj_standard(:,i),Npix,Np_full), w, 1, 0, [Npix,Npix], 2);
    reco_NUFFT_standard(:,:,i) =  FT' * (rawdata_standard(:,:,i));
end

% now reconstruct the off resonant scan
FT = NUFFT(reshape(k_traj_standard_off,Npix,Np_full), w, 1, 0, [Npix,Npix], 2);
reco_NUFFT_standard_off = FT' * rawdata_standard_off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Normalize data             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Normalize the data...');

CEST_standard = abs(reco_NUFFT_standard) ./ ...
    abs(repmat(reco_NUFFT_standard_off, [1 1 Nimages]));

CEST_UCI = abs(reco_S_on_sub) ./ abs(reco_S_off_sub);

% plot normalized data
subplot(1,2,1)
ShowImages(CEST_standard,[],[],[0 1]);
title('normalized standard CEST')

subplot(1,2,2)
ShowImages(CEST_UCI,[],[],[0 1]);
title('normalized UCI')

display('These are the normalized data sets as obtained by dividing pixel-wise' );
display('with an off-resonant scan.');
display('***Press button to continue***'); waitforbuttonpress;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Read in ROIs       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load filenames of ROIs
NRois = 3;
for i = 1:NRois
    RoiPath{i} = ['ROI', num2str(i), '.roi'];    
end
  
% interate over ROIs ...
for k = 1:NRois
    
    % read in current ROI
    file = fopen(RoiPath{k}, 'r');
    
    % read in first entry
    tempx=(fgetl(file)); 
    tempy=(fgetl(file));
    xpos(1,k) = str2num(tempx);
    ypos(1,k) = str2num(tempy);
    count = 1;
    
    % read in remaining entries
    while 1 
        count = count + 1;
        tempx=(fgetl(file)); 
        tempy=(fgetl(file));
        if tempx == -1; break; end      % break if EOF
        xpos(count,k) = str2num(tempx);
        ypos(count,k) = str2num(tempy);
    end
    
    % get ROI mask for current ROI   
    RoiMask(:,:,k) = ...
        roipoly(CEST_UCI(:,:,1), xpos(xpos(:,k) ~= 0, k),ypos(ypos(:,k) ~= 0, k));       
        
    % iterate over images ...
    for i = 1:Nimages
        
        % get UCI spectrum
        currentImage = abs(CEST_UCI(:,:,i));
        [mean_UCI(i,k), sdev_UCI(i,k)] = ...
            EvalROI(currentImage, RoiMask(:,:,k)); 

        % get standard CEST spectrum
        currentImage=abs(CEST_standard(:,:,i));
        [mean_standard(i,k), sdev_std(i,k)] = ...
            EvalROI(currentImage, RoiMask(:,:,k)); 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Plot CEST spectra   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalize spectra to mean intensity of points at 33.7 ppm (see Fig. 4)
divisor = mean(mean_UCI(8,:)) / mean(mean_standard(8,:));
mean_UCI = mean_UCI / divisor;

% plot CEST spectra
subplot(1, 3, 1:2);
Legend = {};
count = 1;
for k=1:NRois
        hold all
        p(k) = plot(ragefreqs, mean_UCI(:,k), 'o');
        hold off
        Legend{count} = 'UCI';
        count=count + 1;
end
for k=1:NRois
        hold on
        plot(reffreqs, mean_standard(:,k),'x:', 'color', ...
            get(p(k), 'color'), 'MarkerSize', 9);
        hold off
        Legend{count} = 'standard';
        count = count + 1;
end

legend(Legend, 'Location', 'SouthEast')
xlabel('chemical shift [ppm]')
ylabel('normalized signal')
title('CEST spectra for B_1 = 15 \mu T')
% show ROIs
subplot(1, 3, 3)
imagesc(abs(reco_NUFFT_standard_off(:,:,1)));colormap gray; 
axis image; axis off

hold on
for k=1:NRois
    DrawROI(xpos(xpos(:,k)~=0,k),ypos(ypos(:,k)~=0,k),get(p(k),'color'),2);
end
hold off
title('regions of interest (ROIs)')

display('The UCI spectra agree well with the spectra obtained with standard CEST.');
