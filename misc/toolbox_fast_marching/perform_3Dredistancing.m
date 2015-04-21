function D1 = perform_3Dredistancing(D, options)

% perform_redistancing - redistance a function
%
%   D1 = perform_redistancing(D, options);
%
%   Compute a signed distance function D1 that have the same 0 level set as D.
%   You can turn off interpolation by using options.use_interpolation=0.
%
%   Note that the distance function is computed with 1 pixel = distance of
%   1, so the overall image range over a 1...n=size(D,1)
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
use_interpolation = getoptions(options, 'use_interpolation', 1);

[height, width,depth] = size(D);

n = max([height,width,depth]);
tmpD = ones(n,n,n)*min(D(:));
tmpD(round((n-height)/2+1):round((n+height)/2),round((n-width)/2+1):round((n+width)/2),round((n-depth)/2+1):round((n+depth)/2)) =D;
D = tmpD;


% horizontal
P1 = D(1:end-1,:,:); P2 = D(2:end,:,:);
P = (P1.*P2)<0;
d = abs(P1-P2); d(d<eps) = 1;
v1 = abs(P1)./d;
v2 = abs(P2)./d;
Ah1 = zeros(n,n,n);
Ah1(1:n-1,:,:)=P;
Ah2 = zeros(n,n,n);
Ah2(2:n,:,:) = P ;
Ah = (Ah1 + Ah2)>0;
Vh1 = zeros(n,n,n);
Vh1(1:n-1,:,:) =v1;
Vh2 = zeros(n,n,n);
Vh2(2:n,:,:)=v2;
Vh = max(Vh1, Vh2);
% vertical
P1 = D(:,1:end-1,:); P2 = D(:,2:end,:);
P = (P1.*P2)<0;
d = abs(P1-P2); d(d<eps) = 1;
v1 = abs(P1)./d;
v2 = abs(P2)./d;
Av1 = zeros(n,n,n);
Av1(:,1:n-1,:)=P;
Av2 = zeros(n,n,n);
Av2(:,2:n,:) = P ;
Av = (Av1 + Av2)>0;
Vv1 = zeros(n,n,n);
Vv1(:,1:n-1,:) =v1;
Vv2 = zeros(n,n,n);
Vv2(:,2:n,:)=v2;
Vv = max(Vv1, Vv2);


% depth
P1 = D(:,:,1:end-1); P2 = D(:,:,2:end);
P = (P1.*P2)<0;
d = abs(P1-P2); d(d<eps) = 1;
v1 = abs(P1)./d;
v2 = abs(P2)./d;
Ad1 = zeros(n,n,n);
Ad1(:,:,1:n-1)=P;
Ad2 = zeros(n,n,n);
Ad2(:,:,2:n) = P ;
Ad = (Ad1 + Ad2)>0;
Vd1 = zeros(n,n,n);
Vd1(:,:,1:n-1) =v1;
Vd2 = zeros(n,n,n);
Vd2(:,:,2:n)=v2;
Vd = max(Vd1, Vd2);

V = zeros(n,n,n);
I = find(Ah>0);
V(I) = Vh(I);
I = find(Av>0);
V(I) = max(V(I),Vv(I));
I = find(Ad>0);
V(I) = max(V(I),Vd(I));

I = find(V~=0);
[x,y,z] = ind2sub(size(D),I); 
start_points = [x(:)'; y(:)';z(:)'];
vals = V(I);

% with interpolation
options.nb_iter_max = Inf;
if use_interpolation
    options.values = vals/n;
else
    options.values = [];    
end
D1 = perform_fast_marching(ones(n,n,n), start_points, options);
D1 = D1*n;
D1(D<0) = -D1(D<0);
D1 = D1(round((n-height)/2+1):round((n+height)/2),round((n-width)/2+1):round((n+width)/2),round((n-depth)/2+1):round((n+depth)/2)) ;