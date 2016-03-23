function imin=imgMarkCrossHair(imin,c1c2,fmtstr)
% Cross hair is 2 lines of value 1, intersecting at c1c2.
% Input img can be 3D
if nargin<3,fmtstr='solid'; end
switch fmtstr
    case 'dashed'
        imin(c1c2(1),1:2:end)=1; % mark with dashed lines
        imin(1:2:end,c1c2(2))=1; % mark with dashed lines
    otherwise
        imin(c1c2(1),:)=1; % mark with solid lines
        imin(:,c1c2(2))=1; % mark with solid lines
end