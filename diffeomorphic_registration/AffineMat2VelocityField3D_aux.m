function [log_def_x, log_def_y, log_def_z] = AffineMat2VelocityField3D_aux(AffMat3D, GridSize)
%function [log_def_x, log_def_y, log_def_z] = AffineMat2VelocityField3D_aux(AffMat3D, GridSize)

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

t = AffMat3D(1:3,4);
M = AffMat3D(1:3,1:3);
c = transpose(GridSize/2);

t = t - M*c + c;

T = [M, t; [0 0 0 1]];
logT = logm(T);

L = logT(1:3,1:3);
v = logT(1:3,4);

[X, Y, Z] = ndgrid(1:GridSize(1), 1:GridSize(2), 1:GridSize(3));

temp_stack = [X(:)'; Y(:)';Z(:)'];
log_def_x = L(1,:)*temp_stack + v(1);
log_def_y = L(2,:)*temp_stack + v(2);
log_def_z = L(3,:)*temp_stack + v(3);

log_def_x = real(reshape(log_def_x, GridSize));
log_def_y = real(reshape(log_def_y, GridSize));
log_def_z = real(reshape(log_def_z, GridSize));

