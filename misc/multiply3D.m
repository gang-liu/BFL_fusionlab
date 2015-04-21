% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
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

function mat3D = multiply3D(input_matrix1, input_matrix2)

% input is 3 x 3 x N
% output is a 3 x 3 x N matrix where each entry is the product of each 3 by
% 3 of the entries

mat3D = cat(1,sum(squeeze(input_matrix1(1,:,:)).*squeeze(input_matrix2(:,1,:))), ...
    sum(squeeze(input_matrix1(2,:,:)).*squeeze(input_matrix2(:,1,:))), ...
    sum(squeeze(input_matrix1(3,:,:)).*squeeze(input_matrix2(:,1,:))), ...
    sum(squeeze(input_matrix1(1,:,:)).*squeeze(input_matrix2(:,2,:))), ...
    sum(squeeze(input_matrix1(2,:,:)).*squeeze(input_matrix2(:,2,:))), ...
    sum(squeeze(input_matrix1(3,:,:)).*squeeze(input_matrix2(:,2,:))), ...
    sum(squeeze(input_matrix1(1,:,:)).*squeeze(input_matrix2(:,3,:))), ...
    sum(squeeze(input_matrix1(2,:,:)).*squeeze(input_matrix2(:,3,:))), ...
    sum(squeeze(input_matrix1(3,:,:)).*squeeze(input_matrix2(:,3,:))));

mat3D = reshape(mat3D, [3 3 size(mat3D,2)]);
return;
    