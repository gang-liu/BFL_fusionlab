function AffMat3D = AffineParams2Mat_aux(params)
%function AffMat3D = AffineParams2Mat_aux(params)

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

logsx = params(1);
logsy = params(2);
logsz = params(3);
rx = params(4);
ry = params(5);
rz = params(6);
tx = params(7);
ty = params(8);
tz = params(9);

S = diag([exp(logsx), exp(logsy), exp(logsz)]);
    
    Rx = [1       0        0; ...
          0 cos(rx) -sin(rx); ...
       0 sin(rx) cos(rx)];
   
    Ry = [cos(ry) 0 sin(ry); ...
          0       1       0; ...
          -sin(ry) 0  cos(ry)];
      
    Rz = [cos(rz) -sin(rz) 0; ...
        sin(rz) cos(rz) 0; ...
        0       0       1];
t = [tx; ty; tz];

R = Rx*Ry*Rz;
AffMat3D = [S *R, t];
AffMat3D = [AffMat3D; 0 0 0 1];
