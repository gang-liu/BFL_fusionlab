%BFL_labelfusion_DEMO.m

SBJ_CELL = {'IBSR_01_mri_norm.mgz', 'IBSR_02_mri_norm.mgz', 'IBSR_03_mri_norm.mgz', 'IBSR_04_mri_norm.mgz'};
DATA_DIR = '../data/';
Output_dir = '../temp/';

% this makes the registration algorithm output the progress of the
% algorithm -- ie the cost function value at each iteration 
% by default this would be zero - so no progress would be printed while
% performing the computationally expensive registration
options.verbose = 1;

options.num_multires = 4;
options.numiter_cell = cell(options.num_multires,1);
options.reg_weight = 20;

for i = 3:options.num_multires
    options.numiter_cell{i,1} = 150;
end
options.numiter_cell{2,1} = 30;
options.numiter_cell{1,1} = 5;
options.rigidFlag = 1; %%% This flag (if nonzero) tells the registration algorithm to pre-register w/ an affine transformation model


%%% this parameter determines the "nonlinearity/flexibility" of the registration 
%typically it should be an integer greater than zero
%the smaller it is, the more nonlinear, deformable the registration model
%is
options.sigma_diff = 2;

%the following parameter determines the weighting in the weighted local
%label fusion -- if it is very high (eg approaches infinitiy) then all
%training subjects get the same weight. smaller values make the training
%images that have a local similar appearance to the test image get weighed
%more.
options.sigma_labelfusion = 5;

for ii= 1:length(SBJ_CELL)

    if (isfield(options, 'TRAINING_IMAGE_CELL'))
        options = rmfield(options, {'TRAINING_IMAGE_CELL', 'TRAINING_SEG_CELL', 'rigidParamsCell'});
    end
    
    %% let's identify training subjects that are the "closest" to the test
    %% subjects
    [RELEV_TRAIN_CELL, RELEV_TRAIN_SEG_CELL, affine_params_cell] = BFL_identify_relevant_training_aux(SBJ_CELL{ii}, DATA_DIR, options);

    options.TRAINING_IMAGE_CELL = RELEV_TRAIN_CELL;
    options.TRAINING_SEG_CELL = RELEV_TRAIN_SEG_CELL;
    options.rigidParamsCell = affine_params_cell;
    BFL_labelfusion_3D({SBJ_CELL{ii}}, DATA_DIR, Output_dir, options)
    
end
