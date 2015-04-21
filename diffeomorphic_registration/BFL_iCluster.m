%function [cluster_indx] = BFL_iCluster(SBJ_CELL, DATA_DIR, Output_dir, params)
function [cluster_indx] = BFL_iCluster(SBJ_CELL, DATA_DIR, Output_dir, params)

% input:
% *SBJ_CELL : a cell with filenames of subject image data
% * DATA_DIR: the directory where the subject image data resides
% * Output_dir: the directory where the results will be written. 
%   There are four sets of results: 
%      (i) warp files, which are .mat files with names like [Output_dir '/' SBJ_CELL{1,i} '.' CurTemplateWarpName '.mat']
%      (ii) warped subjects: these are warped subject images (resampled on
%      to the template grid) -- the names are as follows: [OUTPUT_DIR '/'
%      SBJ_CELL{1,i} Sbj_filename_postfix_out]
%      (iii) the cluster template(s): For each cluster we save the cluster mean (template) as [Output_dir '/TemplateSbj_Clust' num2str(cluster_idx) '.nii']
%       (iv) cluster_indx.txt -- a file with subject name and cluster_indx
%       in each row
% *params: structure with fields
%   mask_vol:  the mask volume in which image-to-image distance will be
%   computed -- at voxels with nonzero mask_vol value
%   reg_weight : <scalar> this effectively determines the step size in the
%   optimization (gradient descent). a higher value will yield smaller
%   steps. (default: 50)
%   initialWarpFileName : the filenames (.mat)for the initial warps. Default will assume [];
%   preregistered_flag: nonzero indicates that the images have been
%   pre-registered - so there is no need to initially groupwise register
%   all images
%   num_outer_iter : number of iteretions of the outside loop (the loop
%   that iterates between computing the "template image" (average subject)
%   and registering all subjects with the template
%   num_inner_iter : number of iterations/pyramid scale/registration
%   initial_sigma_diff: initial sigma for the warp regularization (see
%   sigma_diff option in pairwise_reg3D.m)
%   final_sigma_diff : final sigma for the warp regularization (see
%   sigma_diff option in pairwise_reg3D.m)
%   num_multires : number of resolutions in the multires pyramid
%   Sbj_filename_postfix_out = e.g.
%   '-dwi-filt-Ed-nhdr_masked_Baseline-rgd.nii.gz' or '.nii.gz' -- This
%   will be added to the name of the subject when saving the
%   warped/resampled subjects



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


if (~isfield(params, 'preregistered_flag'))
    preregistered_flag = 0;
else
    preregistered_flag = params.preregistered_flag;
end


if (~isfield(params, 'initialWarpFileName'))
    initialWarpFileName = [];
else
    initialWarpFileName = params.initialWarpFileName;
end


if (preregistered_flag ~= 0 && isempty(initialWarpFileName))
    
    
    

