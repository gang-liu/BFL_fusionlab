%function BFL_groupwise_reg3D(SBJ_CELL, DATA_DIR, Output_dir, params)
function BFL_groupwise_reg3D(SBJ_CELL, DATA_DIR, Output_dir, params)

% input:
% *SBJ_CELL : a cell with filenames of subject image data
% * DATA_DIR: the directory where the subject image data resides
% * Output_dir: the directory where the results will be written. 
%   There are three sets of results: 
%      (i) warp files, which are .mat files with names like [Output_dir '/' SBJ_CELL{1,i} '.' CurTemplateWarpName '.mat']
%      (ii) warped subjects: these are warped subject images (resampled on
%      to the template grid) -- the names are as follows: [OUTPUT_DIR '/'
%      SBJ_CELL{1,i} Sbj_filename_postfix_out]
%      (iii) the template(s): after each iteration, we save the templates as
%      .mat files ([Output_dir '/Template' num2str(iter) '.mat']). The very
%      last template will be saved as [Output_dir '/TemplateSbj.nii']
% *params: structure with fields
%   reg_weight : <scalar> this effectively determines the step size in the
%   optimization (gradient descent). a higher value will yield smaller
%   steps. (default: 50) if this is not set properly (which should be done
%   for each image type since it depends on the intensity range in the
%   images), then the optimization may converge prematurely.
%   num_outer_iter : number of iteretions of the outside loop (the loop
%   that iterates between computing the "template image" (average subject)
%   and registering all subjects with the template
%   num_inner_iter : number of iterations/pyramid scale/registration
%   initial_sigma_diff: initial sigma for the warp regularization (see
%   sigma_diff option in pairwise_reg3D.m)
%   final_sigma_diff : final sigma for the warp regularization (see
%   sigma_diff option in pairwise_reg3D.m) 
%   if this is infinity then only do
%   affine registration
%   
%   num_multires : number of resolutions in the multires pyramid
%   Sbj_filename_postfix_out = e.g.
%   '-dwi-filt-Ed-nhdr_masked_Baseline-rgd.nii.gz' or '.nii.gz' -- This
%   will be added to the name of the subject when saving the
%   warped/resampled subjects
%   initial_template_sbj_name = name of subject who will serve as initial
%   template guess


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


NumSbj = length(SBJ_CELL);

if (~exist(Output_dir, 'dir'))
    mkdir(Output_dir);
end

if (~isfield(params, 'reg_weight'))
    
    reg_weight = 50;
else
    reg_weight = params.reg_weight;
end

if (~isfield(params, 'initial_template_sbj_name'))
    
    initial_template_sbj_name = [];
else
    initial_template_sbj_name = params.initial_template_sbj_name;
end

if (~isfield(params, 'verbose'))
    
    params.verbose = 0;
end

if (~isfield(params, 'renormalize_warps_flag'))
    
    renormalize_warps_flag = 1;
else
    renormalize_warps_flag = params.renormalize_warps_flag;
    
end




if (~isfield(params, 'num_outer_iter'))
    
    num_outer_iter = 6;
else
    num_outer_iter = params.num_outer_iter;
end

if (~isfield(params, 'num_inner_iter'))
    num_inner_iter = 200;
else
    num_inner_iter = params.num_inner_iter;
end


if (~isfield(params, 'sigma_diff'))
    
    final_sigma_diff = 2;
else
    final_sigma_diff = params.sigma_diff;
end

if (~isfield(params, 'num_multires'))
    
    num_multires = 4;
else
    num_multires = params.num_multires;
end


if (~isfield(params, 'Sbj_filename_postfix_out'))
    
    Sbj_filename_postfix_out = '_warped2template.nii';
else
    Sbj_filename_postfix_out = params.Sbj_filename_postfix_out;
end


options.invcon_flag = 1;
options.use_jacobian = 1;
options.fw_weight = 0;
options.gradient_type= 0;
options.sigma_fluid =  0;
options.max_step =  10;

options.numiter =  num_inner_iter;
options.num_multires = num_multires;
options.reg_weight = reg_weight;

if (~isfield(params, 'numiter_cell'))
    options.numiter_cell = cell(options.num_multires, 1);

    for  i = 1:options.num_multires
        options.numiter_cell{i,1} = options.numiter;
    end

    options.numiter_cell{2,1} = 10;
    options.numiter_cell{1,1} = 2;
end
options.verbose = params.verbose;
range = (initial_sigma_diff - final_sigma_diff)/(num_outer_iter - 1);

final_sigma_diff = initial_sigma_diff;

if (range > 0)

    sigma_diff_vec = initial_sigma_diff:-range:final_sigma_diff;
else
    sigma_diff_vec = initial_sigma_diff*ones(num_outer_iter, 1);
end
CurTemplateWarpName = 'CurTemplateWarp';
%compute initial template
if (params.verbose)
    display('Computing initial template...');
end

if (isempty(initial_template_sbj_name))
    template_vol = compute_template_aux(SBJ_CELL, DATA_DIR, Output_dir);

    if (params.verbose)
        display('done...');
    end
else
    vol1 = MRIread([DATA_DIR '/' initial_template_sbj_name]);
    template_vol = vol1.vol;
end

save([Output_dir '/Template0.mat'], 'template_vol');
for iter = 1:num_outer_iter
    if (params.verbose)
        display(['Pairwise registrations... iteration: ' num2str(iter)]);
    end
    options.sigma_diff = sigma_diff_vec(iter);
    
    options.num_multires = max(5 - iter, 2);
    
    if (iter < num_outer_iter)
        options.min_level = 2;
    else
        options.min_level = 1;
    end
    
    load([Output_dir '/Template' num2str(iter-1) '.mat']);
    for i = 1:NumSbj
        if (iter ~= 1)
            load([Output_dir '/' SBJ_CELL{1,i} '.' CurTemplateWarpName '.mat']);
            options.log_def_x = log_def_x;
            options.log_def_y = log_def_y;
            options.log_def_z = log_def_z;
        end
        if (params.verbose)
            display(['Registering subj : ' SBJ_CELL{1,i}]);
        end
        
        if (iter > 1)
            options.rigidFlag = 0;
        end
        vol1 = MRIread([DATA_DIR '/' SBJ_CELL{1,i}]);
        if ~isinf(final_sigma_diff)
            [log_def_x, log_def_y, log_def_z, stats, wm] = ...
                BFL_pairwise_reg3D((vol1.vol), (template_vol), options);
        else
            [affineMat] = ...
                    BFL_pairwise_affine_reg3D((vol1.vol), (template_vol), options);
            [log_def_x, log_def_y, log_def_z] = AffineMat2VelocityField3D_aux(affineMat, size(vol1.vol));
        end
        
        save([Output_dir '/' SBJ_CELL{1,i} '.' CurTemplateWarpName '.mat'], ...
            'log_def_x', 'log_def_y', 'log_def_z');
        clear log_def_x log_def_y log_def_z stats wm
    end
    
    if (renormalize_warps_flag)
        if (params.verbose)
            display('Renormalizing warps...');
        end
        renormalize_warps(SBJ_CELL, Output_dir, CurTemplateWarpName);
    end
    if (params.verbose)
            display('Updating template...');
    end
    template_vol = compute_template_aux(SBJ_CELL, DATA_DIR, Output_dir, CurTemplateWarpName);
    save([Output_dir '/Template' num2str(iter) '.mat'], 'template_vol');

end
OUTPUT_DIR = Output_dir;

i = 1;
output_template = MRIread([DATA_DIR '/' SBJ_CELL{1,i} ]);
output_template.vol = template_vol;
MRIwrite(output_template, [OUTPUT_DIR '/TemplateSbj.nii']);

for i = 1:NumSbj
    load([Output_dir '/' SBJ_CELL{1,i} '.' CurTemplateWarpName '.mat']);
   
    vol1 = MRIread([DATA_DIR '/' SBJ_CELL{1,i}]);
    [def_x, def_y, def_z] = velocityfieldexp(-log_def_x, -log_def_y, -log_def_z);
    vol1.vol = warpimage(single(vol1.vol), def_x, def_y,def_z);
    vol1.vol(isnan(vol1.vol)) = 0;
    temp_cnt = sum([SBJ_CELL{1,i} Sbj_filename_postfix_out] == '/');
    
    if (temp_cnt == 0)
        MRIwrite(vol1, [OUTPUT_DIR '/' SBJ_CELL{1,i} Sbj_filename_postfix_out]);
    else
        jj = find([SBJ_CELL{1,i} Sbj_filename_postfix_out] == '/');
        tmp_str = [SBJ_CELL{1,i} Sbj_filename_postfix_out];
        for ii = jj
            if (~exist([OUTPUT_DIR '/' tmp_str(1:ii-1)], 'dir'))
                disp(['Making directory: ' OUTPUT_DIR '/' tmp_str(1:ii-1)]);
                mkdir([OUTPUT_DIR '/' tmp_str(1:ii-1)]);
            end
        end
        MRIwrite(vol1, [OUTPUT_DIR '/' SBJ_CELL{1,i} Sbj_filename_postfix_out]);
    end
        %need to check whether directory exists
end

