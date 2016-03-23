function [xy, peaks, acum]=findLVctrCine(img,radii,varargin)
% 
% [xy, peaks, acum]=findLVctrCine(img,radii,'npeaks',3,'movname',movname)
% 
% INPUT
% 
%     im: 3D image of size [n1 n2 n3], each page is a frame.
% 
%     radii: e.g. 15:25, all radius values to check, in pxs.
% 
%     'nhoodxy', 'nhoodxy': Num of pixels to smooth the accumulator along xy and radius. Default=3.
%     
%     'npeaks': Top N circles to search for per frame. Defaults to 1.
% 
%     'avgAcum':  If FALSE [default]: run peak detection for each time
%         frame of accumulator. If TRUE: run peak detection after accumulator
%         has been averaged over time frames.
% 
%     'movname': a debug cine is saved in this name. Omit to save no cines.
% 
% OUTPUT
% 
%     xy: [xind; yind], in px, final verdict of LV center, from upper left corner, i.e. dim2,dim1 indeces.
%     
%     peaks: [3, nradii, nframe]: Along col: x, y, radius. Along row: N top circles. Along page: cine frames.
% 
%     acum: [n1, n2, nradius, n3]: Hough accumulator values
%     
%     

% test: 
%     [xy, peaks, acum]=findLVctrCine(im3d,10:25,'npeaks',1,'avgAcum',1,'movname','001 sl006');
%     figure; imagescn(acum,[],[],[],4);    
% 
%     mean(acum,4);
%     figure; imagescn(mean(acum,4),[],[],[],4);
% 

n1=size(img,1);
n2=size(img,2);
n3=size(img,3);
nr=length(radii);

[npeaks, varargin]=findStripArg(varargin,'npeaks',1);
[avgAcum, varargin]=findStripArg(varargin,'avgAcum',false);
[nhoodxy, varargin]=findStripArg(varargin,'nhoodxy',3);
[nhoodr, varargin]=findStripArg(varargin,'nhoodr',3);
[movname, varargin]=findStripArg(varargin,'movname','');

difIm=mostDiffImg(img); %%%%
% difIm=abs(img-repmat(mean(img,3),[1 1 n3]));

acum=zeros(n1,n2,nr,n3);
if ~avgAcum, peaks=zeros(3,npeaks,n3); end % every col: xind, yind, radius. [px]
saveMov=~isempty(movname);

% debug outputs: , img1, img2, img3
if saveMov
    img1=zeros(size(img)); % post-threhold
    img2=img1; % post erosion
end

fprintf('%s: running frame: ',mfilename);
for i3=1:n3
    fprintf('%d ',i3);
    
    imIn=LG_normalize(abs(difIm(:,:,i3)));
    imIn=im2bw(imIn,graythresh(imIn));    
    
%     imIn=LG_normalize(img(:,:,i3));
    
%     imHough=edge(imIn,'canny');    
    imHough=imerode(imIn,strel('disk',1));
%     imHough=imerode(imIn,strel('square',3));
%     imHough=edge(imHough,'canny');    
%     imHough=edge(imHough,'sobel'); 
    
%     Manual gradient doesn't work as well as Canny edge
%     [gx, gy]=gradient(imIn);
%     imIn=LG_normalize(gx.^2+gy.^2);
%     imHough=thresh_img(imIn,prctile(imIn(:),90)); %im2bw(imIn,graythresh(imIn));
    
    acum(:,:,:,i3) = circle_hough(imHough, radii, 'same', 'normalise');  %   sum(h_out,3);
    if ~avgAcum
        peaks(:,:,i3) = circle_houghpeaks(acum(:,:,:,i3), radii, 'nhoodxy', nhoodxy, 'nhoodr', nhoodr, 'npeaks', npeaks);
    end
    if saveMov
        img1(:,:,i3)=imIn;
        img2(:,:,i3)=imHough;
        % img3(:,:,i3)=double(imHough);
    end
end
fprintf('Done\n');

% debug:
%     pks=circle_houghpeaks(mean(acum,4), radii, 'nhoodxy', 3, 'nhoodr', 3, 'npeaks', 3);

if avgAcum
    peaks=circle_houghpeaks(mean(acum,4), radii, 'nhoodxy', nhoodxy, 'nhoodr', nhoodr, 'npeaks', npeaks);
    xy=round(median(peaks(1:2,:),2)); % returns [x;y]
else
    % xy=round(mean(mean(peaks(1:2,:,:),2),3)); % returns [x;y]
    xy=round(median(median(peaks(1:2,:,:),2),3)); % returns [x;y]
end

% saves a 4-up movie
if saveMov
    %%
    img_normd=LG_scale(img);
    img_normd2=img_normd;
    
    % mark circle on each frame:
    for i3=1:n3
        if avgAcum % need this IF branch because dim3 of peaks will have diff sizes.
            for ir=1:npeaks
                [x, y] = circlepoints(peaks(3,ir));
                ilin=sub2ind([n1 n2 n3],y+peaks(2,ir),x+peaks(1,ir),ones(size(x)).*i3);
                img_normd2(ilin)=1;
            end    
        else
            for ir=1:npeaks
                [x, y] = circlepoints(peaks(3,ir,i3));
                ilin=sub2ind([n1 n2 n3],y+peaks(2,ir,i3),x+peaks(1,ir,i3),ones(size(x)).*i3);
                img_normd2(ilin)=1;
            end    
        end
    end
    
    % mark center for cine as crosshair on all frames:
    rx=2;
    img_normd2(xy(2)-rx:xy(2)+rx,xy(1),:)=1; % mark each line of crosshair separately
    img_normd2(xy(2),xy(1)-rx:xy(1)+rx,:)=1;
    img_4up=[img_normd,double(img1);double(img2),img_normd2];
    if n3>1
        dumpMP4_img3D(img_4up,movname);
    else
        imwrite(LG_normalize(img_4up),[movname,'.png']);
    end
end

%% junk
% Plot onto axis to draw the circle:
%                 imagesc(img2(:,:,i3)); colormap gray; hold on;  axis off; daspect([1 1 1]);
%                 for ir=1:npeaks
%                     [x, y] = circlepoints(peaks(3,ir,i3));
%                     plot(x+peaks(1,ir,i3),y+peaks(2,ir,i3),'r-','linewidth',2);
%                 end

% if plotReport
%     figure_LG('Debug plot',[40 40 660 760]);
%         % one iteration per row:
%         nrowdisp=5;
%         for i3=1:n3
%             scrollsubplot(nrowdisp,3,(i3-1)*3+1);
%                 imagesc(img1(:,:,i3)); colormap gray; axis off; daspect([1 1 1]);
%             scrollsubplot(nrowdisp,3,(i3-1)*3+2);
%                 imagesc(img2(:,:,i3)); colormap gray; axis off; daspect([1 1 1]);
%             scrollsubplot(nrowdisp,3,(i3-1)*3+3);
%                 for ir=1:npeaks
%                     [x, y] = circlepoints(peaks(3,ir,i3));
%                     ilin=sub2ind([n1 n2 n3],y+peaks(2,ir,i3),x+peaks(1,ir,i3),ones(size(x)).*i3);
%                     img3(ilin)=0.5; %img3(y+peaks(2,ir,i3),x+peaks(1,ir,i3),i3)=0.5;
%                 end
%                 imagesc(img3(:,:,i3)); colormap gray; hold on;  axis off; daspect([1 1 1]);
%         end
% end


