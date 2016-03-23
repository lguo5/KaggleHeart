function score=maskRegBldInMyo(msk_bld,msk_myo,varargin)
% Scores each discontiguous blood region by how much myocardium is
% surrounding it, and fills the region with its score.
% 
% INPUT
% 
%     msk_bld,msk_myo: 2D masks, should have same size. The two masks should not intercept.
% 
%     'normalize': (defaults to FALSE) if TRUE, will divide region score by pixel count of region.
% 
% OUTPUT
% 
%     score: large is "more surrounded", will be scaled to 0~1. (Or remains all 0 if no score was accumulated.)
% 
% 
[normalizeScore,varargin]=findStripArg(varargin,'normalize',false);

n1=size(msk_bld,1);
n2=size(msk_bld,2);
score=zeros(n1,n2);

L=bwlabel(msk_bld);
nreg=max(L(:));

% if input msk_bld is all FALSE, then the loop will just be skipped.
for ireg=1:nreg
    msk=L==ireg;
    distmap=bwdist(msk);
    score(msk)=sum(distmap(msk_myo));
    if normalizeScore,
        score(msk)=score(msk)./sum(msk(:));
    end
end

% scale region scores to 0.1 and (1.1 - smallest score)
if any(score(:)~=0),
    score=LG_scale(score);
end
score=1.1-score;
score(~msk_bld)=0; % nuke none-blood bgd

% test:
%         maskRegBldInMyo(msk3d_bld(:,:,1),msk3d_myo(:,:,1),'normalize',false);
