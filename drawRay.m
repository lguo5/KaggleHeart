function acum=drawRay(acum,d1,d2,deg,intens,maxR,bidir)
% Draws multiple rays (additive overlay) in an accumulator matrix.
% 
% INPUTS
% 
%     acum: accumulator matrix to work on, usually an all-0 matrix need not be.
%     
%     d1, d2: Nx1 each, [d1, d2] is ray origin, in dim1 and dim2 subscripts. This point is NOT drawn to avoid repeated pileup.
% 
%     deg: Nx1: angles of rays in degrees, -180~180 CCW from the positive x-axis (dim2)
%
%     intens: Nx1: intensities of rays.
% 
%     maxR: num of pixels: radius of each ray drawn, e.g. 30. Use [] to draw up to image border.
% 
%     bidir: logical scalar, whether to draw each ray in both directions (making a traversing line).
% 
% 
% OUTPUT
% 
%     acum: accumulator after drawing.
%     
%     
if nargin<7 || isempty(bidir), bidir=false; end
if nargin<6, maxR=[]; end
    
n1=size(acum,1);
n2=size(acum,2);

d1=d1(:);
d2=d2(:);
deg=deg(:);
intens=intens(:);

% if bidir
%     deg=[deg;180+deg];
%     d1=[d1;d1];
%     d2=[d2;d2];
%     intens=[intens;intens];
% end
nray=size(d1,1);

maxlen=ceil(sqrt(n1^2+n2^2)); % distance in num of pxs. worst case distance
if isempty(maxR),    maxR=maxlen;
else    maxR=min(maxlen,maxR);
end

if bidir
    lind=-maxR:maxR;
else
    lind=0:maxR; % indeces along line
end
lind=lind(:);

for ind=1:nray
    % mark pixels in dim 1 and dim 2 to be affected
    x1=d1(ind)+round(lind.*  sind(deg(ind)) );
    x2=d2(ind)+round(lind.*(-cosd(deg(ind))));
    keepmask=~((x1>n1)|(x1<1)|(x2>n2)|(x2<1));
    x1=x1(keepmask);
    x2=x2(keepmask);
%     icut=find((x1>n1)|(x1<1)|(x2>n2)|(x2<1),1,'first');
%     if ~isempty(icut)
%         x1=x1(1:icut-1);
%         x2=x2(1:icut-1);
%     end
    iOn=sub2ind([n1 n2],x1,x2);
    iOn=unique(iOn); % needed to prevent value piling on a pixel
    acum(iOn)=acum(iOn)+intens(ind);
end

% test:
%     drawRay(zeros(300,400),repmat(30,[4 1]),repmat(30,[4 1]),[-150;-40;50;120],1:4,100);    
%     drawRay(zeros(300,400),repmat(30,[4 1]),repmat(30,[4 1]),[-150;-40;50;120],1:4,100,true);      

















