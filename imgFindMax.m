function [d1,d2,d3,maxval]=imgFindMax(im)
% Finds dim1 dim2 dim3 positions of the max intensity in a 2D or 3D image. 
%     [d1,d2,d3,maxval]=matrixFindMax(im2d)
% Note that d3 will be 1 if input is a 2D image.
% 
% 

% n1=size(im,1);
% n2=size(im,2);
% n3=size(im,3);

[maxval,maxind]=max(im(:));
[d1,d2,d3]=ind2sub(size(im),maxind);
