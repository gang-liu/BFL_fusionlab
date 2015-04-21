function [affineMat, affineParams] = ...
    BFL_pairwise_affine_reg3D_multires(fix_im, mov_im, options)
%function [affineMat, affineParams] = ...
%    BFL_pairwise_affine_reg3D_multires(fix_im, mov_im, options)

% Uses Gauss-Newton optimization
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
    
    options.num_of_used_voxels = 5000;
end

if ~(isfield(options, 'foreground_mask_im'))
   options.foreground_mask_im  = fix_im;
end

if ~isfield(options, 'MaxIter')
    options.MaxIter = 200;
end

if ~isfield(options, 'ParamsTol')
    options.ParamsTol = 0.01;
end

if ~isfield(options, 'num_multires')
    numOfLevels = 4;
else
    numOfLevels = options.num_multires;
end

pyramid1 = cell(numOfLevels, 1);

pyramid2 = cell(numOfLevels, 1);

pyramid1{1,1} = fix_im;

pyramid2{1,1} = mov_im;



for level = 2:1:numOfLevels



    size1 = size(pyramid1{level-1,1});

    size1 = size1 - mod(size1,2);

    pyramid1{level,1} = pyramid1{level-1,1}(1:2:size1(1),1:2:size1(2), 1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),1:2:size1(2),1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(1:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(1:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(1:2:size1(1),2:2:size1(2),2:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),2:2:size1(2),2:2:size1(3));



    pyramid1{level,1} = pyramid1{level,1}/8;



    size2 = size(pyramid2{level-1,1});

    size2 = size2 - mod(size2,2);

    pyramid2{level,1} = pyramid2{level-1,1}(1:2:size2(1),1:2:size2(2), 1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),1:2:size2(2),1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(1:2:size2(1),2:2:size2(2),1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(1:2:size2(1),1:2:size2(2),2:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(1:2:size2(1),2:2:size2(2),2:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),1:2:size2(2),2:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),2:2:size2(2),1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),2:2:size2(2),2:2:size2(3));



    pyramid2{level,1} = pyramid2{level,1}/8;



end

options.initial_affineParams = [0 0 0 0 0 0 0 0 0]';

for level = numOfLevels: -1: 1

    options.foreground_mask_im = pyramid1{level,1};
    %[affineMat, affineParams] = ...
    %    BFL_pairwise_affine_reg3D_fast(pyramid1{level,1}, pyramid2{level,1}, options);

    [affineMat, affineParams] = ...
        BFL_pairwise_affine_reg3D(pyramid1{level,1}, pyramid2{level,1}, options);
    options.initial_affineParams = affineParams;
    options.initial_affineParams(end-3:end) = affineParams(end-3:end)*2;
    
end

affineMat = AffineParams2Mat_aux(affineParams);

