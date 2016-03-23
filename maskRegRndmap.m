function [rMap, info]=maskRegRndmap(msk,varargin)
% Finds roundness metric for each discontiguous region in a 2D binary
% image, then fills that region with its roundness as intensity.
% 
% e.g. [rMap,info]=maskRegRndmap(msk,'convex',true,'normalize',true)
% 
% INPUT 
% 
%     msk: 2D binary image of discrete regions.
%     
%     'convex', true/false. If true, finds roundness of the conv hull of each region, instead of roundness of region itself.
% 
%     'normalize', true/false. whether to normalize roundness to num of pixels when painting a region with its roundness
% 
% 
% OUTPUT
% 
%     rMap: accumulator of roundness, same size as msk, larger value is higher roundness.
%     
%     info: nObjx1 struct array
% 

[normalizeRndMap,varargin]=findStripArg(varargin,'normalize',false);
[useConvex,varargin]=findStripArg(varargin,'convex',false);

n1=size(msk,1);
n2=size(msk,2);
[B,L] = bwboundaries(msk,'noholes'); % Every cell element of B is the coordinates of boundary points: Nx2: [dim1 dim2]. There will be nObj elements (number of detected objects) in B.

stats = regionprops(L,{'Area','Centroid','ConvexArea','ConvexHull'});

rMap=zeros(size(msk)); % accumulator; higher is rounder

nObj=length(B);
for iObj=1:nObj
    
    if useConvex % convex metrics:
        vtx=stats(iObj).ConvexHull;
        tmp.mask=poly2mask(vtx(:,1),vtx(:,2),n1,n2);    
        tmp.area=stats(iObj).ConvexArea;
        tmp.circumf=sum(sqrt(sum([diff(vtx,1,1);vtx(1,:)-vtx(end,:)].^2,2)));
    else
        tmp.mask=L==iObj;% poly2mask(B{iObj}(:,2),B{iObj}(:,1),n1,n2); % notice B's [dim1 dim2] is [y x]
        tmp.area=stats(iObj).Area;
        tmp.circumf=sum(sqrt(sum(diff(B{iObj},1,1).^2,2))); % connect boundary points w lines then find line lenghts
    end
    tmp.centroid=stats(iObj).Centroid; % !! centroid in convex mode is still the original region's
    tmp.roundness=4*pi*tmp.area/tmp.circumf^2;
    
    if normalizeRndMap
        incr=tmp.roundness/sum(L(:)==iObj);
    else
        incr=tmp.roundness;
    end
    rMap(L==iObj)=rMap(L==iObj)+incr; % !! even in convex mode, roundness is added to the original shape of the region
    
    if iObj==1
        info=repmat(tmp,[nObj 1]);
    end
    info(iObj)=tmp;
end


%% junk

% before making convex a function input option:
%     [B,L] = bwboundaries(msk,'noholes'); 
%     stats = regionprops(L,{'Area','Centroid','ConvexArea','ConvexHull'});
%     tmp.roundnessMap=zeros(size(msk)); % accumulator; higher is rounder
%     tmp.cvxRoundnessMap=zeros(size(msk)); % accumulator; higher is rounder
% 
%     nObj=length(B);
%     for iObj=1:nObj
% 
%         tmp.circumf = sum(sqrt(sum(diff(B{iObj}).^2,2))); % connect boundary points w lines then find line lenghts
%         tmp.mask=L==iObj;% poly2mask(B{iObj}(:,2),B{iObj}(:,1),n1,n2); % notice B's [dim1 dim2] is [y x]
%         tmp.area=stats(iObj).Area;
%         tmp.centroid=stats(iObj).Centroid;
%         tmp.roundness=4*pi*tmp.area/tmp.circumf^2;
%         if normalizeRndMap
%             incr=tmp.roundness/sum(L(:)==iObj);
%         else
%             incr=tmp.roundness;
%         end
%         tmp.roundnessMap(L==iObj)=tmp.roundnessMap(L==iObj)+incr;
% 
%         % convex metrics:
%         tmp.cvxArea=stats(iObj).ConvexArea;
%         vtx=stats(iObj).ConvexHull;
%         tmp.cvxCircumf=sum(sqrt(sum([diff(vtx);vtx(end,:)-vtx(1,:)].^2,2)));
%         tmp.cvxRoundness=4*pi*tmp.cvxArea/tmp.cvxCircumf^2;
%         tmp.cvxHullMask=poly2mask(vtx(:,1),vtx(:,2),n1,n2);    
%         if normalizeRndMap
%             incr=tmp.cvxRoundness/sum(L(:)==iObj); % normalize to original region, NOT its conv hull
%         else
%             incr=tmp.cvxRoundness;
%         end
%         tmp.cvxRoundnessMap(L==iObj)=tmp.cvxRoundnessMap(L==iObj)+incr;
% 
%         if iObj==1
%             out=repmat(tmp,[nObj 1]);
%         end
%         out(iObj)=tmp;
%     end
% 
