%function template_vol = compute_template_aux(SBJ_CELL, DATA_DIR, Output_dir, CurTemplateWarpName)
function template_vol = compute_template_aux(SBJ_CELL, DATA_DIR, Output_dir, CurTemplateWarpName)

%%% THIS IS AN AUXILARY FUNCTION AND IS NOT MEANT TO BE CALLED DIRECTLY.
%%% THE FUNCTION TO BE CALLED IS: groupwise_reg3D

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



use_jac_flag = 0;
NumSbj = length(SBJ_CELL);

for i = 1:NumSbj
   
    vol1 = MRIread([DATA_DIR '/' SBJ_CELL{1,i}]);
    
    if(i == 1)
        template_vol = zeros(size(vol1.vol), 'single');
    end
    
    if (nargin >= 4)
        
        load([Output_dir '/' SBJ_CELL{1,i} '.' CurTemplateWarpName '.mat']);
        [def_x, def_y, def_z] = velocityfieldexp(-log_def_x, -log_def_y, -log_def_z);
        vol1.vol = warpimage(single(vol1.vol), single(def_x), single(def_y),single(def_z));
        vol1.vol(isnan(vol1.vol)) = 0;
        
        if (use_jac_flag)
           weight_vol = deffieldjacobiandeterminant(def_x, def_y, def_z);
        end
    end
    
    if (~use_jac_flag)
        weight_vol = ones(size(vol1.vol), 'single');
    end
    
    if (i == 1)
        total_weight_vol = zeros(size(vol1.vol), 'single');
    end
        
    template_vol = weight_vol.*vol1.vol + template_vol;
    total_weight_vol = weight_vol + total_weight_vol;
    
end

template_vol = template_vol./total_weight_vol;