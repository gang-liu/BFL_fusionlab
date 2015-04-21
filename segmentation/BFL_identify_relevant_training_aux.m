%function [RELEV_TRAIN_CELL, RELEV_TRAIN_SEG_CELL, affine_params_cell] = BFL_identify_relevant_training_aux(SBJ, DATA_DIR, options)
function [RELEV_TRAIN_CELL, RELEV_TRAIN_SEG_CELL, affine_params_cell] =  BFL_identify_relevant_training_aux(SBJ, DATA_DIR, options)

% This is a small module that will identify the relevant training subjects
% by performing a quick affine registration between the test subject and
% the training images -- and picking the training images that are "closest"
% to the test image
% input:
% *SBJ: filename of subject image data
% * DATA_DIR: the directory where the subject image data resides
% * Output_dir: the directory where the results will be written.
% * options: a structure with following potential fields:
%
%   *TRAINING_DATA_DIR: the directory where the training data data resides
%   *TRAINING_IMAGE_CELL: filenames for training images (eg)
%   *TRAINING_SEG_CELL: filenames for training label images (ie manual
%   segmentations)
%   background_label = <int> the index of the background label (default: 0)

%=========================================================================
%  Brain Fuse Lab 
%  Copyright (c) 2010-2011 Mert R. Sabuncu
%  A.A. Martinos Center for Biomedical Imaging, 
%  Mass. General Hospital, Harvard Medical School
%  CSAIL, Mass. Institute of Technology 
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================

if (nargin<3)
    options = [];
end

if (~isfield(options,'TRAINING_DATA_DIR'))
    options.TRAINING_DATA_DIR = '../data/';
end

if (~isfield(options,'TRAINING_IMAGE_CELL'))
    options.TRAINING_IMAGE_CELL = {'IBSR_01_mri_norm.mgz' 'IBSR_03_mri_norm.mgz' 'IBSR_05_mri_norm.mgz' ...
        'IBSR_07_mri_norm.mgz' 'IBSR_09_mri_norm.mgz' 'IBSR_11_mri_norm.mgz' 'IBSR_13_mri_norm.mgz' ...
        'IBSR_15_mri_norm.mgz' 'IBSR_17_mri_norm.mgz' 'IBSR_02_mri_norm.mgz' ...
        'IBSR_04_mri_norm.mgz' 'IBSR_06_mri_norm.mgz' 'IBSR_08_mri_norm.mgz' ...
        'IBSR_10_mri_norm.mgz' 'IBSR_12_mri_norm.mgz' 'IBSR_14_mri_norm.mgz' ...
        'IBSR_16_mri_norm.mgz' 'IBSR_18_mri_norm.mgz'...
        };
end

if (~isfield(options,'NumRelTrainingSubjects'))
    options.NumRelTrainingSubjects = 8;
end

if (~isfield(options,'TRAINING_SEG_CELL'))
    options.TRAINING_SEG_CELL = {'IBSR_01_man_seg.mgz' 'IBSR_03_man_seg.mgz' 'IBSR_05_man_seg.mgz' ...
        'IBSR_07_man_seg.mgz' 'IBSR_09_man_seg.mgz' 'IBSR_11_man_seg.mgz' 'IBSR_13_man_seg.mgz' ...
        'IBSR_15_man_seg.mgz' 'IBSR_17_man_seg.mgz' 'IBSR_02_man_seg.mgz' ...
        'IBSR_04_man_seg.mgz' 'IBSR_06_man_seg.mgz' 'IBSR_08_man_seg.mgz' ...
        'IBSR_10_man_seg.mgz' 'IBSR_12_man_seg.mgz' 'IBSR_14_man_seg.mgz' ...
        'IBSR_16_man_seg.mgz' 'IBSR_18_man_seg.mgz'...
        };
end


if (~isfield(options,'background_label'))
    options.background_label = 0;
end

if (~isfield(options,'labels'))
    seg2 = MRIread([options.TRAINING_DATA_DIR '/' options.TRAINING_SEG_CELL{1}]);

    labels = unique(seg2.vol(:));
else
    labels = options.labels;
    
end

labels = setdiff(labels, options.background_label);
NumTrainingSubjects = length(options.TRAINING_IMAGE_CELL);

% let's first generate mask

seg2 = MRIread([options.TRAINING_DATA_DIR '/' options.TRAINING_SEG_CELL{1}]);
mask = zeros(size(seg2.vol));

for i = 1:NumTrainingSubjects
        
     seg2 = MRIread([options.TRAINING_DATA_DIR '/' options.TRAINING_SEG_CELL{i}]);
     
     for l = 1:length(labels)
         mask(seg2.vol == labels(l)) = 1;
     end
end

temp_cost_vec = zeros(NumTrainingSubjects,1);
RELEV_TRAIN_CELL = cell(options.NumRelTrainingSubjects,1);
RELEV_TRAIN_SEG_CELL= cell(options.NumRelTrainingSubjects,1);
affine_params_cell = cell(options.NumRelTrainingSubjects,1);
tmp_vol = MRIread([DATA_DIR '/' SBJ]);

fix_vol = tmp_vol.vol;

display('Identifying relevant subset of training subjects...');
affineParams_stack = zeros(NumTrainingSubjects, 9);
for i = 1:NumTrainingSubjects

    %%% you don't want to be including the test subject in the training
    %%% set
    if ~(strcmp(options.TRAINING_IMAGE_CELL{i}, SBJ))
        display(['Processing training subj. ' num2str(i) ' out of ' num2str(NumTrainingSubjects)]);
        tmp_vol = MRIread([options.TRAINING_DATA_DIR '/' options.TRAINING_IMAGE_CELL{i}]);

        mov_vol = tmp_vol.vol;

        tic
        [affineMat_est, affineParams_est] = ...
            BFL_pairwise_affine_reg3D_multires(fix_vol, mov_vol);
        toc
        affineParams_stack(i,:) = affineParams_est;
        warped_mov_vol = resample_im_affine3D_aux(mov_vol, affineParams_est);
        warped_mov_vol(isnan(warped_mov_vol)) = 0;

        temp_cost_vec(i) = mean((fix_vol(mask > 0) - warped_mov_vol(mask > 0)).^2);
    else
        temp_cost_vec(i) = NaN;
    end
end

[Y, tmp_ind] = sort(temp_cost_vec);

for i = 1:options.NumRelTrainingSubjects
    RELEV_TRAIN_CELL{i} = options.TRAINING_IMAGE_CELL{tmp_ind(i)};
    RELEV_TRAIN_SEG_CELL{i} = options.TRAINING_SEG_CELL{tmp_ind(i)};
    affine_params_cell{i} = affineParams_stack(tmp_ind(i),:);
end

return

