function [msk3d_myo,msk3d_bld, intens_myo,intens_bld, msk3d_hrt,im3d_mopx]=cineClusterBldMyo(im3d,varargin)
% INPUT
%     im3d: each page is a cine frame
%     'hrtMskDilateR': e.g. 5 pixel. Dilates heart mask to include more myocardium. Defaults to no dilation.
%     'mskTrimSaveFrxn': defaults to 0.8. Small motion regions are trimmed (to save at least this much total area) before their convex hull is made.
% OUTPUT
%     msk3d_myo,msk3d_bld: 3D masks, same size as im3d, that marks myocardium and blood
%     intens_myo,intens_bld: determined mean intensity of myo and blood
%     msk3d_hrt: 3D masks, same size as im3d, that marks the heart region. All pages are the same!
%     im3d_mopx: 3D image, same size as im3d, that shows the beating heart, in a black background.
% 

[motionRegTrimSaveFrxn,varargin]=findStripArg(varargin,'motionRegTrimSaveFrxn',0.8);
[hrtMskDilateR,varargin]=findStripArg(varargin,'hrtMskDilateR',0); % defaults to no dilation

n1=size(im3d,1);
n2=size(im3d,2);
n3=size(im3d,3);

% motion pixels
imdif=std(im3d,0,3); % LG_scale(std(im3d,0,3)); 
msk_mopx=im2bw(imdif,graythresh(imdif));
msk_mopx=medfilt2(msk_mopx,[5 5]);  %%%%
msk_mopx=imfill(msk_mopx,'holes'); % to add mass to closed regions, usually LV

% Deleting smallest regions, save top 80 or 90% area, then make a convex hull of what's left:
regions_mopx=bwconncomp(msk_mopx,4); %%%%
regions_mopx_npx=cellfun(@(x)length(x),regions_mopx.PixelIdxList);
[regions_mopx_npx_sorted,sortInd]=sort(regions_mopx_npx,'descend');
regions_mopx_npx_sorted_cumsum=cumsum(regions_mopx_npx_sorted)./sum(regions_mopx_npx_sorted);
iEnough=find(regions_mopx_npx_sorted_cumsum>motionRegTrimSaveFrxn,1,'first'); %%%%
msk_mopx=false(size(msk_mopx));
msk_mopx(cell2mat(transpose(regions_mopx.PixelIdxList(sortInd(1:iEnough)))))=true;
msk_hrt=bwconvhull(msk_mopx);
if hrtMskDilateR>0
    msk_hrt=imdilate(msk_hrt,strel('disk',hrtMskDilateR));
end

% makes a "heart-only" 3D image
msk3d_hrt=repmat(msk_hrt,[1 1 n3]);
im3d_mopx=zeros(size(im3d));
im3d_mopx(msk3d_hrt)=im3d(msk3d_hrt);

% k-means clustering of intensity
[classNum,C]=kmeans(im3d(msk3d_hrt),2);
im_label=zeros(size(im3d));
im_label(msk3d_hrt)=classNum;
[C,sortInds]=sort(C);
msk3d_myo=im_label==sortInds(1); % the darker of the two classes, assumed to be myocardium
msk3d_bld=im_label==sortInds(2);  %  diff(msk3d_bld,1,3);    std(msk3d_bld,0,3);    

intens_myo=C(1);
intens_bld=C(2);