function [diffImgs,diffMat]=mostDiffImg(img)
% For each frame in a 3D image (a cine), returns its difference image
% with another frame that results in the largest-difference.

n1=size(img,1);
n2=size(img,2);
n3=size(img,3);
% diffImgs=zeros(size(img),'like',img);

img=reshape(img,[n1*n2,n3]);

% calc upper triangular
diffMat=zeros(n3,n3);
for ii=1:(n3-1)
    for jj=(ii+1):n3
        diffMat(ii,jj)=norm(img(:,ii)-img(:,jj));
    end
end
% fold and make full matrix, find max along each row (might as well be cols)
diffMat=diffMat+diffMat';
[~, maxinds]=max(diffMat,[],2);

diffImgs=img-img(:,maxinds);
diffImgs=reshape(diffImgs,[n1,n2,n3]); % abs(diffImgs)
