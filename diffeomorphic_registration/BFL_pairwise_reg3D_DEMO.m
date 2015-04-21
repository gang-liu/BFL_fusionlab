%BFL_pairwise_reg3D_DEMO.m

%SBJ_CELL = {'IBSR_01_mri_norm.mgz', 'IBSR_10_mri_norm.mgz'};
%DATA_DIR = '../data/';

%vol1 = MRIread([DATA_DIR '/' SBJ_CELL{1}]);
%vol2 = MRIread([DATA_DIR '/' SBJ_CELL{2}]);


clear

vol_orignal01 = MRIread('02singlelabel_1.nii.gz');
vol_orignal02 = MRIread('12singlelabel_1.nii.gz');
fix_vol_affine = vol_orignal01.vol; % only translate
mov_vol_affine = vol_orignal02.vol;
[affineMat, affineParams] = ...
    BFL_pairwise_affine_reg3D_multires(fix_vol_affine, mov_vol_affine);

vol1 = MRIread('02singlelabel_1.nii.gz');
vol2 = MRIread('12singlelabel_1.nii.gz');
[log_def_x, log_def_y, log_def_z] = AffineMat2VelocityField3D_aux(affineMat, size(fix_vol_affine));
[def_x, def_y, def_z] = velocityfieldexp(log_def_x, log_def_y, log_def_z);
warped_vol2 = warplabelimage(vol2.vol, def_x, def_y, def_z);
warped_vol2(isnan(warped_vol2))=0;
vol2.vol=warped_vol2;



fix_vol = vol1.vol;
mov_vol = vol2.vol;
options.verbose = 1;
options.num_multires = 4;
options.rigidFlag = 0; % No pre-affine transformation
options.reg_weight = 1;


options.numiter_cell = cell(options.num_multires,1);
for i = 1:options.num_multires
    options.numiter_cell{i} = 100;
end
options.numiter_cell{2} = 100;
options.numiter_cell{1} = 100;

options.sigma_diff = 2;
[log_def_x, log_def_y, log_def_z, stats, warped_mov_im, backwarped_fix_im] = ...
    BFL_pairwise_reg3D(fix_vol, mov_vol, options, affineMat);

checkerboard_imagesc(squeeze(fix_vol(:,:,128)),squeeze(mov_vol(:,:,128)), 1)
figure(1), title('Fixed Image and Moving Image');
checkerboard_imagesc(squeeze(fix_vol(:,:,128)),squeeze(warped_mov_im(:,:,128)), 2)
figure(2), title('Fixed Image and Warped Moving Image');



