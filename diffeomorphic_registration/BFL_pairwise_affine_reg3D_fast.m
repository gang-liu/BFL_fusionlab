function [affineMat, affineParams] = ...
    BFL_pairwise_affine_reg3D_fast(fix_im, mov_im, options)
%function [affineMat, affineParams] = ...
%    BFL_pairwise_affine_reg3D_fast(fix_im, mov_im, options)

% Uses Gauss-Newton optimization
%input: 
% fix_im: a 3D double matrix (e.g. MRI volume of size 256x256x256)
% mov_im: a 3D double matrix (e.g. MRI volume of size 256x256x256)
%
% NOTE: This algorithm assumes fix_im and mov_im are of the same size!!!
GridSize = size(fix_im);
warning off all
if nargin<3
    options = [];
end

if ~isfield(options,'num_of_used_voxels')
    
    num_of_used_voxels = 5000;
else
    num_of_used_voxels = options.num_of_used_voxels;
end

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
   
    initial_affineParams = [0 0 0 0 0 0 0 0 0]';
    
else
    initial_affineParams = options.initial_affineParams;
end

[X, Y, Z] = ndgrid(1:GridSize(1), 1:GridSize(2), 1:GridSize(3));

pnts = [X(voxel_flat_ind_list_used) Y(voxel_flat_ind_list_used) Z(voxel_flat_ind_list_used)];

fix_im_vals = fix_im(voxel_flat_ind_list_used);
mov_im_gradX = zeros(size(mov_im));
mov_im_gradY = zeros(size(mov_im));
mov_im_gradZ = zeros(size(mov_im));

n=3;
xvals=(-n:n);
gausskernel = exp(-xvals.^2/2);
gausskernel = gausskernel/sum(gausskernel);


mov_im_sm = imfilter(mov_im, gausskernel, 'replicate');
mov_im_sm = imfilter(mov_im_sm, gausskernel', 'replicate');
mov_im_sm = imfilter(mov_im_sm, shiftdim(gausskernel,-1), 'replicate');

mov_im_gradX(2:end,:,:)= diff(mov_im_sm, 1,1);
mov_im_gradY(:,2:end,:)= diff(mov_im_sm, 1,2);
mov_im_gradZ(:,:,2:end)= diff(mov_im_sm, 1,3);

mov_im_gradX_vals = mov_im_gradX(voxel_flat_ind_list_used);
mov_im_gradY_vals = mov_im_gradY(voxel_flat_ind_list_used);
mov_im_gradZ_vals = mov_im_gradZ(voxel_flat_ind_list_used);

%params_limit = [0.3 0.3 0.3 pi/4 pi/4 pi/4 size(fix_im,1)/4 size(fix_im,2)/4 size(fix_im,3)/4];
%[fout,gout,jac] = @(params)cost_fun(params, fix_im_vals, mov_im, mov_im_gradX, mov_im_gradY, mov_im_gradZ, pnts);
[affineParams,histout,costdata] = gaussn(initial_affineParams, ...
    @(params)cost_fun(params, fix_im_vals, mov_im, mov_im_gradX_vals, mov_im_gradY_vals, mov_im_gradZ_vals, pnts),...
    ParamsTol,MaxIter);

affineMat = AffineParams2Mat_aux(affineParams);

return;


function [fout,gout,jac]=cost_fun(params, fix_im_vals, mov_im, mov_im_gradX_vals, mov_im_gradY_vals, mov_im_gradZ_vals, pnts)

AffPtx = AffineWarp3D_aux(pnts,params,size(mov_im));

warped_mov_im_vals = interpn(mov_im, AffPtx(:,1), AffPtx(:,2), AffPtx(:,3), 'linear', 0);

r = (fix_im_vals - warped_mov_im_vals);

GradMat = AffineGrad3D_aux(pnts,params,size(mov_im));

jac = zeros(length(r), 9);

for i =1:9
   jac(:,i) = -sum([mov_im_gradX_vals mov_im_gradY_vals mov_im_gradZ_vals].*squeeze(GradMat(:,:,i)), 2);    
end
fout = r'*r/2;
gout = jac'*r;

return;



