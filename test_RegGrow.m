
%% load
[im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_6'); % most basal   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_7'); % basal   
% [im3d_raw, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_12');    


n1=size(im3d_raw,1);
n2=size(im3d_raw,2);

im3d=im3d_raw(round(n1/4)+1:round(n1/4*3),round(n2/4)+1:round(n2/4*3),:);
n1=size(im3d,1);
n2=size(im3d,2);
n3=size(im3d,3);

%% find blood, myo intens in motion column (mocol)

% I=im3d(:,:,1); % I(70,60)  
% ctr_d1=70;  % case 1, sax_12, frame 1
% ctr_d2=60;  % case 1, sax_12, frame 1
% ctr_rad=30; % case 1, sax_12, frame 1

I=im3d(:,:,11); % case 1, sax_6
ctr_d1=65;
ctr_d2=58;
ctr_rad=30;

% How to generate motion mask -----------

% threshold diff image:
% msk=im2bw(LG_normalize(abs(im3d(:,:,10)-im3d(:,:,1))),.5);

% square FOV:
msk=false(size(im3d(:,:,1)));
msk(ctr_d1-ctr_rad:ctr_d1+ctr_rad,ctr_d2-ctr_rad:ctr_d2+ctr_rad)=true;
% ---------

msk_3d=repmat(msk,[1 1 n3]);
mocol_img=im3d;
mocol_img(~msk_3d)=0;

mocol_val=im3d(msk_3d); % col vec.   figure; hist(mocol_val,50);   
val_myo=prctile(mocol_val,30);
val_bld=prctile(mocol_val,90);
val_mid=(val_myo+val_bld)/2;
val_gap=val_bld-val_myo;

%%
[mask_reg, mask_nei, rec_intens, rec_mean, rec_traj]=...
    regiongrowing_LG(I,ctr_d1,ctr_d2,'pxval',[val_myo+val_gap/4 inf]);

figure;
imshowpair(I,mask_reg,'blend','scaling','independent');

%% 
I_msk=zeros(size(I),'like',I);
I_msk(mask_reg)=I(mask_reg);
I_msk_2=sigmoid((I_msk-val_gap)./(val_gap/2).*4.6); % sig(-4.6)=0.01  sig(4.6)=0.99
fprintf('num of effective LV px: %d\n',round(sum(I_msk_2(:))));

%%
% J = regiongrowing(im3d(:,:,1),70,60,95); 