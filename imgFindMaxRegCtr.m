function [d1d2,imgMarked,mskMaxReg,stats]=imgFindMaxRegCtr(img,intensFrxn)
% For a 2D image where the max value may be found in (one or more) region:
% >> Finds an intensity thresh by some fraction of intensity dynamic range (NOT
% population), then masks the regions above it.
% >> Returns the centroid of the largest of such regions: d1d2 for [dim1,dim2],
% and also the mask (same size of img) of that region.
% 
% INPUTS
%     img: 2D image
%     intensFrxn: fraction (default 90) of intensity dynamic range
% OUTPUTS
%     d1d2: [dim1 dim2] inds of centroid of the max region.
%     imgMarked: input image with crosshairs centered on d1d2.
%     mskMaxReg: mask of the max region, same image size as input image.
%     
n1=size(img,1);
n2=size(img,2);

intensMin=min(img(:));
intensRange=max(img(:))-intensMin;
mskRaw=img>intensMin+intensFrxn*intensRange;

L=bwlabel(mskRaw);
stats=regionprops(L,{'Area','Centroid'});

[~,indMaxArea]=max(cell2mat({stats.Area}));

mskMaxReg=L==indMaxArea;
d1d2=flip(round(stats(indMaxArea).Centroid)); % centroid output is x,y, need flip

d1d2(d1d2<1)=1;
if d1d2(1)>n1, d1d2(1)=n1; end;
if d1d2(2)>n2, d1d2(2)=n2; end;

imgMarked=img;
imgMarked(d1d2(1),:)=1;
imgMarked(:,d1d2(2))=1;

