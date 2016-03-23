function [msk,msknpx,d1f,d2f]=findLVarea_RegGro(im2d,d1,d2,intensBld,intensBndBld,boxHW, varargin)
% INPUTS
%     d1,d2: initial landing point: dim 1, dim2.
%     intensBld: Expected central intens of blood
%     intensBndBld: [lower higher] bounds of acceptible blood intensity, not necessarily centered around intensBld, as intensBndBld maybe asymmetric.
%     varargin: passed onto the region growth routine.
% OUTPUTS
%     d1f,d2f: final landing point: dim 1, dim2.
% 
% If initial landing point is within acceptible blood intensity bounds,
% start region grow. Otherwise search within a box (side=2*boxHW px) for a
% landing point whose intens is most similar to blood. If intensity at this
% point is still outside acceptible intensity bounds, return all-false
% mask.
% 
n1=size(im2d,1);
n2=size(im2d,2);

% if im2d(d1,d2)>=intensBld-intensTolBld  &&  im2d(d1,d2)<=intensBld+intensTolBld
if im2d(d1,d2)>=intensBndBld(1)  &&  im2d(d1,d2)<=intensBndBld(2)
    d1f=d1;
    d2f=d2;
    msk=regiongrowing_LG(im2d,d1f,d2f,varargin{:});
else
    % make wrapped bounding box 
    d1a=max(d1-boxHW,1);
    d1b=min(d1+boxHW,n1);
    d2a=max(d2-boxHW,1);
    d2b=min(d2+boxHW,n2);
    % find min spot:
    diffimgSource=abs(im2d-intensBld);
    diffimg=inf(n1,n2);
    diffimg(d1a:d1b,d2a:d2b)=diffimgSource(d1a:d1b,d2a:d2b);
    [~,imin]=min(diffimg(:));
    [d1f,d2f]=ind2sub([n1 n2],imin);
    % if im2d(d1f,d2f)>=intensBld-intensTolBld  &&  im2d(d1f,d2f)<=intensBld+intensTolBld
    if im2d(d1f,d2f)>=intensBndBld(1)  &&  im2d(d1f,d2f)<=intensBndBld(2)
        msk=regiongrowing_LG(im2d,d1f,d2f,varargin{:});
    else
        msk=false(n1,n2);
    end
end
msknpx=sum(msk(:));
