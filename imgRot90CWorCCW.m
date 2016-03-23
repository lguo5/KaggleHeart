function [imout, isCW]=imgRot90CWorCCW(imA, imB)
% Rotates imA 90deg CW or CCW (decided automatically) to best match imB.
% imA and B can be 3D.
% 
imA_orig=imA;
imA=LG_scale(imA);
imB=LG_scale(imB);
imACW= flip(permute(imA,[2 1 3]),2);
imACCW=flip(permute(imA,[2 1 3]),1);
prodCW= imACW.*imB;
prodCCW=imACCW.*imB;

isCW=sum(prodCW(:))>sum(prodCCW(:));
if isCW,
    imout=flip(permute(imA_orig,[2 1 3]),2);
else
    imout=flip(permute(imA_orig,[2 1 3]),1);
end
