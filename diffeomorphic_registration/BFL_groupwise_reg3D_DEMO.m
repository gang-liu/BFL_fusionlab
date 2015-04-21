%BFL_groupwise_reg3D_DEMO.m

%%% these are the images to be co-registered
SBJ_CELL = {'IBSR_01_mri_norm.mgz', 'IBSR_02_mri_norm.mgz', 'IBSR_03_mri_norm.mgz', ...
    'IBSR_04_mri_norm.mgz', 'IBSR_05_mri_norm.mgz'};


% this is where the image data resides on the disc
DATA_DIR = '../data/';


% this is where we'd like to save our results (the warps as *mat files, and the
% templates -- by default the algorithm saves the intermediate templates
Output_dir = '../temp/';

% this makes the registration algorithm output the progress of the
% algorithm -- ie the cost function value at each iteration 
% by default this would be zero - so no progress would be printed while
% performing the computationally expensive registration
params.verbose = 1;

% this is the image to be selected as the initial template -- note that
% since the warps are renormalized after the completion of all pairwise
% registrations (such that the average deformation is zero across all
% images), the selection of this image is not *too* critical. In other
% words, you shouldn't be too worried about selecting an image of average
% size, etc. 
params.initial_template_sbj_name = 'IBSR_05_mri_norm.mgz' ;

%%% this parameter determines the "nonlinearity/flexibility" of the registration 
%typically it should be an integer greater than zero
%the smaller it is, the more nonlinear, deformable the registration model
%is 
params.sigma_diff = 2;

BFL_groupwise_reg3D(SBJ_CELL, DATA_DIR, Output_dir, params);

for i = 1:5
    load([DATA_DIR '/Template' num2str(i) '.mat']);
    figure, imagesc(squeeze(template_vol(:,:,128))), colormap gray, title(['Template - iter ' num2str(i)]);
end