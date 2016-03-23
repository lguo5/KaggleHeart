function d1d2=maskCentroid(msk)
% Finds the centroid of a binary image, without regard to the connectedness
% among regions.
% 
%     d1d2: [d1 d2] the centroid,  in units of pixels, in [y x] convention
% 
n1=size(msk,1);
n2=size(msk,2);
mass=sum(msk(:));
[g2,g1]=meshgrid(1:n2,1:n1); % meshgrid(1:x,1:y)

m1=msk.*g1;
m2=msk.*g2;

d1d2=[sum(m1(:))/mass,sum(m2(:))/mass];

