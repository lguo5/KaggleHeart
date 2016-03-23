
%% load
[im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_6'); % most basal   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_7'); % basal   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_12'); % mid   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_16'); % apical
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_18'); % most apical
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_38'); % most apical
n1=size(im3d_raw,1);
n2=size(im3d_raw,2);

im3d=im3d_raw(round(n1/4)+1:round(n1/4*3),round(n2/4)+1:round(n2/4*3),:);
n1=size(im3d,1);
n2=size(im3d,2);
n3=size(im3d,3);

%% 
imIn=im3d(:,:,12);
% imIn=imageStacks_echo_1(:,:,20);

msk=false(size(imIn));
msk(65,58)=true;
% msk(56,76)=true; % 

msk=imdilate(msk,strel('disk',4));
figure; imshowpair(imIn,msk,'blend');

%% 
msk2=activecontour(imIn,msk,'Chan-Vese','SmoothFactor',0.8,'ContractionBias',-0.10);
figure; imshowpair(imIn,msk2,'blend');

