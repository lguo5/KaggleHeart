function imout=imgCropOrExpand(imin,n1t,n2t)
% Crops or expands an image to a new dimension [n1t,n2t]. Center of
% original image will be aligned to center of new image. Expanded areas are
% filled with 0s.
% e.g.
%      imout=imgCropOrExpand(rand(200,320),220,340);
% 
% t: target, used on output image
% s: source, used on input image

n1s=size(imin,1);
n2s=size(imin,2);
if islogical(imin)
    imout=false(n1t,n2t);
else
    imout=zeros(n1t,n2t,'like',imin);
end

% i1sa: For source img dim 1, starting index
% i2tb: For target img dim 2, ending index
[i1ta,i1tb,i1sa,i1sb]=calcinds(n1s,n1t);
[i2ta,i2tb,i2sa,i2sb]=calcinds(n2s,n2t);
imout(i1ta:i1tb,i2ta:i2tb)=imin(i1sa:i1sb,i2sa:i2sb);

end

% test:
%     load clown;    X;    
%     imgCropOrExpand(X,100,160);     
%     imgCropOrExpand(X,250,320);     
%     imgCropOrExpand(X,200,370);     
%     imgCropOrExpand(X,250,370);     

function [ita,itb,isa,isb]=calcinds(ns,nt)
% ita: index, for target img, start
% isb: index, for source img, end
    if nt>ns
        ita=max(round((nt-ns)/2),1); % e.g. when ()/2=0.4, change to 1.
        itb=ita+ns-1;
        isa=1;
        isb=ns;
    else
        ita=1;
        itb=nt;
        isa=max(round((ns-nt)/2),1);
        isb=isa+nt-1;
    end

end



%% junk
% if n1t>n1s
%     i1ta=max(round((n1t-n1s)/2),1); % e.g. when ()/2=0.4, change to 1.
%     i1tb=i1ta+n1s-1;
%     i1sa=1;
%     i1sb=n1s;
% else
%     i1ta=1;
%     i1tb=n1t;
%     i1sa=max(round((n1s-n1t)/2),1);
%     i1sb=i1sa+n1t-1;
% end

