function [affineMat, affineParams] = ...
    BFL_pairwise_affine_reg3D(fix_im, mov_im, options)
%function [affineMat, affineParams] = ...
%    BFL_pairwise_affine_reg3D(fix_im, mov_im, options)

%input: 
% fix_im: a 3D double matrix (e.g. MRI volume of size 256x256x256)
% mov_im: a 3D double matrix (e.g. MRI volume of size 256x256x256)
%
% NOTE: This algorithm assumes fix_im and mov_im are of the same size!!!

GridSize = size(fix_im);

if nargin<3
    options = [];
end

if ~isfield(options,'num_of_used_voxels')
    
    num_of_used_voxels = 1000;
else
    num_of_used_voxels = options.num_of_used_voxels;
end


%%% registration will be based on voxels where the foreground_mask_im has
%%% an intensity value greater than zero
if ~(isfield(options, 'foreground_mask_im'))
   options.foreground_mask_im  = fix_im;
end

voxel_flat_ind_list_all = find(options.foreground_mask_im > 0);


if ~isfield(options, 'MaxIter')
    MaxIter = 200;
else
    MaxIter = options.MaxIter;
end

if ~isfield(options, 'ParamsTol')
    ParamsTol = 0.01;
else
    ParamsTol = options.ParamsTol;
end

step_size = round(length(voxel_flat_ind_list_all)/num_of_used_voxels);


voxel_flat_ind_list_used = voxel_flat_ind_list_all(1:step_size:end);

if ~isfield(options, 'initial_affineParams')
   
    initial_affineParams = [0 0 0 0 0 0 0 0 0];
    
else
    initial_affineParams = options.initial_affineParams;
end

[X, Y, Z] = ndgrid(1:GridSize(1), 1:GridSize(2), 1:GridSize(3));

pnts = [X(voxel_flat_ind_list_used) Y(voxel_flat_ind_list_used) Z(voxel_flat_ind_list_used)];

fix_im_vals = fix_im(voxel_flat_ind_list_used);

cost_fun = @(params)ssd_im(params, fix_im_vals, mov_im, pnts);


affineParams = fminsearchbnd(cost_fun, initial_affineParams, [-.3 -.3 -.3 -pi/4 -pi/4 -pi/4 -GridSize(1)/2 -GridSize(2)/2 -GridSize(3)/2], ...
    [.3 .3 .3 pi/4 pi/4 pi/4 GridSize(1)/2 GridSize(2)/2 GridSize(3)/2], optimset('MaxIter', MaxIter, 'TolX', ParamsTol));

affineMat = AffineParams2Mat_aux(affineParams);

return;

function cost = ssd_im(params, fix_im_vals, mov_im, pnts)

AffPtx = AffineWarp3D_aux(pnts,params,size(mov_im));

warped_mov_im_vals = interpn(mov_im, AffPtx(:,1), AffPtx(:,2), AffPtx(:,3), 'linear', 0);

cost = sum((fix_im_vals(:) - warped_mov_im_vals(:)).^2);

return;





