%BFL_pairwise_affine_reg3D_DEMO.m

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
SBJ_CELL = {'IBSR_01_mri_norm.mgz'};
DATA_DIR = '../data/';

vol1 = MRIread([DATA_DIR '/' SBJ_CELL{1}]);

mov_vol = vol1.vol;
fix_vol = resample_im_affine3D_aux(mov_vol, [0 0 0 0 0 0 -20 10 10]); % only translate


tic
[affineMat_est, affineParams_est] = ...
    BFL_pairwise_affine_reg3D_multires(fix_vol, mov_vol);
toc

warped_mov_vol = resample_im_affine3D_aux(mov_vol, affineParams_est);

checkerboard_imagesc(squeeze(fix_vol(:,:,128)),squeeze(mov_vol(:,:,128)), 1)
figure(1), title('Fixed Image and Moving Image');
checkerboard_imagesc(squeeze(fix_vol(:,:,128)),squeeze(warped_mov_vol(:,:,128)), 2)
figure(2), title('Fixed Image and Warped Moving Image');

