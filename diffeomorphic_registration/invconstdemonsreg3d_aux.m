function [log_def_x, log_def_y, log_def_z, stats, warped_mov_im, backwarped_fix_im] = invconstdemonsreg3d_aux(fix_im, mov_im, options)
%function [log_def_x, log_def_y, log_def_z, stats, warped_mov_im, backwarped_fix_im] = invconstdemonsreg3d_aux(fix_im, mov_im, options)

%%% THIS IS AN AUXILARY FUNCTION AND IS NOT MEANT TO BE CALLED DIRECTLY.
%%% THE FUNCTION TO BE CALLED IS: BFL_pairwise_reg3D

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
%      contributors may be used to endorse or promote products derived from
%      this
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



% default options
if ~isfield(options,'numiter')
    options.numiter = 200;
end
if ~isfield(options,'sigma_diff')
    options.sigma_diff = 1;
end
if ~isfield(options, 'reg_weight')
    options.reg_weight = 100;
end
if isfield(options,'log_def_x')
   log_def_x = options.log_def_x;
   log_def_y = options.log_def_y;
   log_def_z = options.log_def_z;
else
   log_def_x = zeros(size(fix_im),class(fix_im));
   log_def_y = zeros(size(fix_im),class(fix_im));
   log_def_z = zeros(size(fix_im),class(fix_im));
end
if ~isfield(options,'verbose')
    options.verbose = 0;
end


if ~isfield(options, 'invcon_flag')
    options.invcon_flag = 1;
end
if ~isfield(options,'sigma_fluid')
    options.sigma_fluid = 0;
end

if options.sigma_diff>0.5
    n=ceil(options.sigma_diff*3);
    x=(-n:n);
    gausskernel_diff = exp(-x.^2/(2*options.sigma_diff^2)) / (options.sigma_diff*sqrt(2*pi));
end

if options.sigma_fluid>0.5
    n=ceil(options.sigma_fluid*3);
    x=(-n:n);
    gausskernel_fluid = exp(-x.^2/(2*options.sigma_fluid^2)) / (options.sigma_fluid*sqrt(2*pi));
else
    gausskernel_fluid = 1;
end

stats.MSE = zeros(options.numiter,1);
stats.backMSE = zeros(options.numiter,1);
stats.harmoEner = zeros(options.numiter,1);
stats.backharmoEner = zeros(options.numiter,1);
stats.negJacRatio = zeros(options.numiter,1);
stats.backnegJacRatio = zeros(options.numiter,1);
fix_im = double(fix_im);
mov_im = double(mov_im);
total_time = 0;

strike = 0;
min_MSE = inf;
for i=1:options.numiter
    try
        if (options.verbose)
            disp(i)
        end
        [log_def_x, log_def_y, log_def_z, stats, up_time] = ...
            make_update(fix_im, mov_im, log_def_x, log_def_y, log_def_z, options, stats, i, gausskernel_fluid, gausskernel_diff);

        total_time = total_time + up_time;
        if (options.verbose)
            disp(stats.MSE(i))
        end
    catch
        err = lasterror;
        disp(err.message)
        break
    end
    
    if ( i > 1)
        if (stats.MSE(i) > (stats.MSE(i - 1) - 1e-8))
            strike = strike + 1;
        else 
            strike = 0;
        end
    end
    
    if (stats.MSE(i) < min_MSE)
        min_MSE = stats.MSE(i);
        log_def_x_min = log_def_x;
        log_def_y_min = log_def_y;
        log_def_z_min = log_def_z;
    end
    
    display(['Min MSE so far: ' num2str(min_MSE) ]);
    
    if (strike > 2)
        
        display('Terminating the optimization...');
        log_def_x = log_def_x_min;
        log_def_y = log_def_y_min;
        log_def_z = log_def_z_min;
        break;
    end 
    
end

stats.total_time = total_time;
if (options.verbose)
    disp(['Elapsed time: ' num2str(total_time) ' seconds']);
end


if ( nargout > 4 )
    [def_x, def_y, def_z] = velocityfieldexp(double(log_def_x), double(log_def_y), double(log_def_z));
    warped_mov_im = warplabelimage(double(mov_im), double(def_x), double(def_y), double(def_z));
    warped_mov_im( isnan(warped_mov_im) ) = 0;
end
if ( nargout > 5 )
    [invdef_x, invdef_y, invdef_z] = velocityfieldexp(-double(log_def_x), -double(log_def_y), -double(log_def_z));
    backwarped_fix_im = warplabelimage(double(fix_im), double(invdef_x), double(invdef_y), double(invdef_z));
    backwarped_fix_im( isnan(backwarped_fix_im) ) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [log_def_x, log_def_y, log_def_z, stats, up_time] = make_update(...
    fix_im, mov_im, log_def_x, log_def_y, log_def_z, options, stats, i, gausskernel_fluid, gausskernel_diff)

tic;

[def_x, def_y, def_z] = velocityfieldexp(log_def_x, log_def_y, log_def_z);
[invdef_x, invdef_y, invdef_z] = velocityfieldexp(-log_def_x, -log_def_y, -log_def_z);
jacdet = deffieldjacobiandeterminant(def_x, def_y, def_z);
inv_jacdet = deffieldjacobiandeterminant(invdef_x, invdef_y, invdef_z);
options.use_jacobian = 1;
if (options.invcon_flag)
    jac_weight = inv_jacdet;
else
    options.use_jacobian = 1;
    jac_weight = 0*inv_jacdet;
end

if (~options.invcon_flag && options.fw_weight)
    [up_x, up_y, up_z] = weightedfwdemonsforces(double(fix_im), double(mov_im), double(def_x), double(def_y), double(def_z), double(invdef_x), double(invdef_y), double(invdef_z), jac_weight, jacdet, options.reg_weight);
else
    [up_x, up_y, up_z] = invcondemonsforces(double(fix_im), double(mov_im), double(def_x), double(def_y), double(def_z), double(invdef_x), double(invdef_y), double(invdef_z), ...
        jac_weight, options.use_jacobian, options.reg_weight);
end

up_time = toc;

% compute stats
warped_mov_im = warplabelimage(mov_im, def_x, def_y, def_z);
idx = isnan(warped_mov_im(:));
stats.MSE(i) = mean( (warped_mov_im(~idx)-fix_im(~idx)).^2 );
stats.harmoEner(i) = deffieldharmonicenergy(def_x, def_y, def_z);

stats.negJacRatio(i) = mean(jacdet(:)<=0);

backwarped_fix_im = warplabelimage(fix_im, invdef_x, invdef_y, invdef_z);
invidx = isnan(backwarped_fix_im(:));
stats.backMSE(i) = mean( (backwarped_fix_im(~invidx)-mov_im(~invidx)).^2 );
stats.backharmoEner(i) = deffieldharmonicenergy(invdef_x, invdef_y, invdef_z);
jacdet = inv_jacdet;
stats.backnegJacRatio(i) = mean(jacdet(:)<=0);

%%%%%


tic;

if 0
if options.sigma_fluid>0.5
    up_x = imfilter(up_x, gausskernel_fluid, 'replicate');
    up_x = imfilter(up_x, gausskernel_fluid', 'replicate');
    up_x = imfilter(up_x, shiftdim(gausskernel_fluid,-1), 'replicate');
    up_y = imfilter(up_y, gausskernel_fluid, 'replicate');
    up_y = imfilter(up_y, gausskernel_fluid', 'replicate');whi
    up_y = imfilter(up_y, shiftdim(gausskernel_fluid,-1), 'replicate');
    up_z = imfilter(up_z, gausskernel_fluid, 'replicate');
    up_z = imfilter(up_z, gausskernel_fluid', 'replicate');
    up_z = imfilter(up_z, shiftdim(gausskernel_fluid,-1), 'replicate');
end
end
log_def_x = log_def_x + up_x;
log_def_y = log_def_y + up_y;
log_def_z = log_def_z + up_z;
if 0
if options.sigma_diff>0.5
    log_def_x = imfilter(log_def_x, gausskernel_diff, 'replicate');
    log_def_x = imfilter(log_def_x, gausskernel_diff', 'replicate');
    log_def_x = imfilter(log_def_x, shiftdim(gausskernel_diff,-1), 'replicate');
    log_def_y = imfilter(log_def_y, gausskernel_diff, 'replicate');
    log_def_y = imfilter(log_def_y, gausskernel_diff', 'replicate');
    log_def_y = imfilter(log_def_y, shiftdim(gausskernel_diff,-1), 'replicate');
    log_def_z = imfilter(log_def_z, gausskernel_diff, 'replicate');
    log_def_z = imfilter(log_def_z, gausskernel_diff', 'replicate');
    log_def_z = imfilter(log_def_z, shiftdim(gausskernel_diff,-1), 'replicate');
end
end
up_time = up_time + toc;