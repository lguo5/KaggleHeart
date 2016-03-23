

%% targets

% Subj 1 -------------------
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_6'); % most basal   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_7'); % basal   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_12'); % mid   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_16'); % apical
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_18'); % most apical
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_38'); % most apical

% Subj 14: Grossly hypertrophied LV --------------
[im3d_raw, dcmInfo]=loadCine('/Users/liheng/Datasets/DSB16/train/14/study/sax_17'); % mid-apical

% crop -------
n1=size(im3d_raw,1);
n2=size(im3d_raw,2);
im3d=im3d_raw(round(n1/4)+1:round(n1/4*3),round(n2/4)+1:round(n2/4*3),:);
% im3d=im3d_raw;

n1=size(im3d,1);
n2=size(im3d,2);
n3=size(im3d,3);

im3d=LG_normalize(im3d,[],[0 99]);


%% time gradient
% imdt=diff(im3d,1,3);  %   abs(imdt)   
[~,~,gz]=gradient(im3d); %    gx.^2+gy.^2+gz.^2;   abs(gz);   
gz2=convn(abs(gz),ones(3,3,3),'same');

acum=zeros(size(im3d));
gmag=zeros(size(im3d));
gdir=zeros(size(im3d));
[a1, a2]=ind2sub([n1,n2],1:(n1*n2));

im4grad=im3d;
% im4grad=im3d-repmat(mean(im3d,3),[1 1 n3]);
% im4grad=gz2;

for ifr=1:n3
    [gmag(:,:,ifr),gdir(:,:,ifr)]=imgradient(im4grad(:,:,ifr));  %  
end
%   gmag.* abs(gz);   
% msk_gmag=gmag>prctile(gmag(:),80);
% gmag=gmag.*gz2;

for ifr=1:n3
    gmag_curr=gmag(:,:,ifr);
    gdir_curr=gdir(:,:,ifr);
    msk=gmag_curr>prctile(gmag_curr(:),80); %     msk=msk_gmag(:,:,ifr);
    % gmag_curr=ones(n1,n2);  % HACK
    acum(:,:,ifr)=drawRay(zeros(n1,n2),a1(msk(:)),a2(msk(:)),gdir_curr(msk(:)),gmag_curr(msk(:)),30,true);
end
%   tmp=mean(acum,3);   % max(tmp(:))

%% kmeans clustering
% [idx,C]=kmeans(im3d(msk_mopx3d),2);
% imtype=zeros(size(im3d));
% imtype(msk_mopx3d)=idx;


%% thresh test
% imIn=LG_normalize(im3d(:,:,10));
% imIn=LG_normalize(abs(imfft(:,:,2)));
% imOut=thresh_img(imIn,.5);   
% imOut=im2bw(imIn,graythresh(imIn)); %   im2bw(imIn,.3)  
% figure; imshow(imOut);

%% candidate diff images

% imdif=mostDiffImg(im3d);
% imdif=diff(im3d,1,3);   %   abs(imdif);   sum(abs(imdif),3);     
% imdif=abs(im3d-repmat(mean(im3d,3),[1 1 n3]));  %   sum(imdif,3)
% imdif=std(im3d,0,3);  % edge(imdif,'canny');  % imgradient(imdif);   

% normim=sum(reshape(im3d,[n1*n2,n3]).^2,1);
% normimdiff=sum(reshape(imdif,[n1*n2,n3]).^2,1);

% imfft=fft(im3d,[],3);  %  abs(imfft)  
% imdiffshft=abs(im3d-circshift(im3d,round(n3/1.5),3));


%% Hough on edges
% 
% % imHough=edge(im3d(:,:,1),'canny',[0.05 .2]);
% % [imHough, th]=edge(abs(imdiffmean(:,:,10)),'canny',[]);
% % [imHough, th]=edge(abs(imfft(:,:,2)),'canny',[0.05 0.2]);
% % edge(abs(imdiffmean(:,:,10)).*imHough,'canny');
% 
% % imIn=LG_normalize(abs(imfft(:,:,2)));
% imIn=LG_normalize(imdif(:,:,10)); % 10: sys. 23: dias
% % imIn=LG_normalize(sum(imdiffmean,3));
% % imHough=edge(im2bw(imIn,graythresh(imIn)),'canny');
% imHough=im2bw(imIn,graythresh(imIn));
% % imHough=imdilate(imHough,strel('square',3));
% % imHough=imerode(imHough,strel('disk',1));
% imHough=imerode(imHough,strel('square',3));
% figure; imshow(imHough);
% 
% radii=15:25;
% h_out = circle_hough(imHough, radii, 'same', 'normalise');  %   sum(h_out,3);    
% peaks = circle_houghpeaks(h_out, radii, 'nhoodxy', 11, 'nhoodr', 5, 'npeaks', 3);
% 
% figure;
% imshow(imHough);
% hold on;
% 
% % plot just 1st circle ----
% % [x, y] = circlepoints(peaks(3,1));
% % plot(x+peaks(1,1), y+peaks(2,1), 'r-','linewidth',2);
% 
% % plot all top circles ---
% for peak = peaks
%     [x, y] = circlepoints(peak(3));
%     plot(x+peak(1), y+peak(2), 'r-','linewidth',2);
% end
% 
% return;


%% Hough on grayscale
% % imHough=im3d(:,:,1);
% % imHough=imdiffmean(:,:,10);
% % imHough=im3d(:,:,10)-im3d(:,:,1);   %    abs(im3d(:,:,10)-im3d(:,:,1));   
% im1=LG_normalize(abs(im3d(:,:,10)-im3d(:,:,1)));
% imHough=double(im2bw(im1,.48)); 
% % imHough=abs(imfft(:,:,2));
% 
% grdthres=(max(imHough(:))-min(imHough(:)))*0.10;
% figure;
% [accum, circen, cirrad, dbg_LMmask] = CircularHough_Grd(imHough, [10 25], grdthres);
%
% th=prctile(accum(:),99);
% figure;
% imshow(accum>th);