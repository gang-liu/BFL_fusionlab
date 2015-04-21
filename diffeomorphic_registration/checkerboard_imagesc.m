function checkerboard_imagesc(im1, im2, fig_id, box_size)

if nargin < 4
    box_size = round(size(im1,1)/10);
end

if (nargin < 3)
    fig_id = figure();
end

figure(fig_id);

composite_im = im1;

box_start_ind1 = 1:box_size:size(im1,1);
box_start_ind2 = 1:box_size:size(im1,2);

box_end_ind1 = box_size-1:box_size:size(im1,1);
box_end_ind2 = box_size-1:box_size:size(im1,2);

if (box_end_ind1(end) ~= size(im1,1))
    
    box_end_ind1 = [box_end_ind1, size(im1,1)];
end
if (box_end_ind2(end) ~= size(im1,2))
    box_end_ind2 = [box_end_ind2, size(im1,1)];
end

for i = 1:2:length(box_start_ind1)
    ii = box_start_ind1(i);
    ii_end = box_end_ind1(i);
    for j = 1:2:length(box_start_ind2)
        jj = box_start_ind2(j);
        jj_end = box_end_ind2(j);
        
        
        composite_im(ii:ii_end, jj:jj_end) = 0.7*im2(ii:ii_end, jj:jj_end);
    end
end

for i = 2:2:length(box_start_ind1)
    ii = box_start_ind1(i);
    ii_end = box_end_ind1(i);
    for j = 2:2:length(box_start_ind2)
        jj = box_start_ind2(j);
        jj_end = box_end_ind2(j);
        
        
        composite_im(ii:ii_end, jj:jj_end) = 0.7*im2(ii:ii_end, jj:jj_end);
    end
end

imagesc(composite_im), colormap gray



