function dt = fast_compute_distance_transform(labIm, label, frameWidth)

if (nargin < 3)
    frameWidth = 5;
end

if (nargin < 2)
    label = 1;
end


I = find(labIm == label);
[X, Y, Z] = ind2sub(size(labIm), I);

Xmin = max(min(X - frameWidth), 1);
Xmax = min(max(X + frameWidth), size(labIm, 1));

Ymin = max(min(Y - frameWidth), 1);
Ymax = min(max(Y + frameWidth), size(labIm, 2));

Zmin = max(min(Z - frameWidth), 1);
Zmax = min(max(Z + frameWidth), size(labIm, 3));

croppedLabIm = (labIm(Xmin:Xmax, Ymin:Ymax, Zmin:Zmax) == label) - 0.5;

fm_options.use_interpolation = 1;
tmp = perform_3Dredistancing(croppedLabIm, fm_options);
dt = min(tmp(:))*ones(size(labIm), 'single');
dt(Xmin:Xmax, Ymin:Ymax, Zmin:Zmax) = tmp;
