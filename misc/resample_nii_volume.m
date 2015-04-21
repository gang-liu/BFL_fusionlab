function new_vol = resample_nii_volume(vol_fixed, vol_mov)

GridSize = vol_fixed.volsize;

[Xind Yind Zind] = ndgrid(1:GridSize(1), 1:GridSize(2), 1:GridSize(3));

%  Creates a structure similar to the FreeSurfer MRI struct
%   defined in mri.h. Times are in ms and angles are in radians.
%   The vox2ras0 matrix is the matrix that converts a 0-based
%   column, row, and slice to XYZ. vox2ras1 is the same with
%   1-based indices. The volume is rows, cols, slices frames,
%   but the vox2ras expects col, row, slice.

R = (vol_fixed.vox2ras1)*[Yind(:) Xind(:) Zind(:) ones(size(Xind(:),1),1)]';

ind2 = transpose((vol_mov.vox2ras1)^(-1)*R);

new_vol = interpn(vol_mov.vol, ind2(:,2), ind2(:,1), ind2(:,3), 'nearest*',0);

new_vol = reshape(new_vol, GridSize);