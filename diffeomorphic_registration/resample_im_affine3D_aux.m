function warped_mov_im = resample_im_affine3D_aux(mov_im, affine_params)
%function warped_mov_im = resample_im_affine3D_aux(mov_im, affine_params)

GridSize = size(mov_im);
[X, Y, Z] = ndgrid(1:GridSize(1), 1:GridSize(2), 1:GridSize(3));

pnts = [X(:) Y(:) Z(:)];
AffPtx = AffineWarp3D_aux(pnts,affine_params,size(mov_im));

warped_mov_im = interpn(mov_im, AffPtx(:,1), AffPtx(:,2), AffPtx(:,3), 'linear', 0);
warped_mov_im = reshape(warped_mov_im, size(mov_im));

