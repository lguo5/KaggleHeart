%% intro

% NEED TO TEST: 
%     SYSTOLIC HEART
%     CASE WHERE ONE SLICE HAS SHIFT

%% load 

% the order of loading and stacking:
% Case 1: has a resp position shift:
path_case='/Users/liheng/Datasets/DSB16/train/1/study';
list_cineFolders={'sax_5';'sax_6';'sax_7';'sax_8';'sax_9';'sax_10';'sax_11';'sax_12';'sax_13';'sax_14';'sax_15';'sax_16';'sax_17';'sax_18'};
cropbox=[85,154,73,152]; % d1a,d1b,d2a,d2b

% % Case 4: 
% path_case='/Users/liheng/Datasets/DSB16/train/4/study';
% list_cineFolders={'sax_5';'sax_6';'sax_7';'sax_8';'sax_9';'sax_10';'sax_11';'sax_12';'sax_13';'sax_14'};
% cropbox=[80,139,50,97]; % d1a,d1b,d2a,d2b

ncine=length(list_cineFolders);

for icine=1:ncine
    [imcurr, dcmInfo]=loadCine([path_case,filesep,list_cineFolders{icine}],'frames',1);
    if icine==1
        n1=size(imcurr,1);
        n2=size(imcurr,2);
        n3=size(imcurr,3);
        im3d_raw=zeros(n1,n2,n3);
    end
    im3d_raw(:,:,icine)=imcurr;
end
im3d_raw=cat(3,im3d_raw(:,:,1),im3d_raw,im3d_raw(:,:,end));
% debug view:    squeeze(im3d(:,73,:));   
% debug view:    permute(im3d,[1 3 2]);        

%% crop
im3d_crop=im3d_raw(cropbox(1):cropbox(2),cropbox(3):cropbox(4),:);
n1=size(im3d_crop,1);
n2=size(im3d_crop,2);
n3=size(im3d_crop,3);
im3d_crop=LG_scale(im3d_crop);

% interp
itpFactor=4;
im3d=interpn(im3d_crop,1:n1,(1:n2)',linspace(1,n3,n3*itpFactor));
% im3d=im3d_crop;
n3=size(im3d,3);
% debug view:
%     squeeze(im3d(:,73,:));   
%     permute(im3d,[1 3 2]);        
%     1-LG_scale(im3d);        

%% thresh
im3d_bw=double((1-im3d)>0.75); % debug view:      permute(im3d_bw,[1 3 2]);        

%% From Snake3D.m:

load testvolume
load SphereMesh  

FV.vertices=FV.vertices./2;

FV.vertices(:,1)=FV.vertices(:,1)+n1/2;
FV.vertices(:,2)=FV.vertices(:,2)+n2/2;
FV.vertices(:,3)=FV.vertices(:,3)+n3/2;

% plot loaded mesh
%   figure; patch(FV,'facecolor',[1 0 0],'facealpha',0.1);
%   figure; plot3(FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3),'s'); camproj perspective

%% run 3D snake
Options=struct;
Options.Verbose=1;
Options.Wedge=2.;
Options.Wline=0.05;
Options.Alpha=0.2;
Options.Beta=0.2;
Options.Kappa=0.9;
Options.Delta=0.1000;
Options.Gamma=0.1000;
Options.Iterations=50;
Options.Sigma1=2;
Options.Sigma2=2;
Options.Lambda=0.8;


% OV=Snake3D(im3d_bw,FV,Options);
OV=Snake3D(im3d,FV,Options);
camproj perspective;  pbaspect([1 1 1]);

%% voxelize mesh
msk_LV=VOXELISE(1:n1, 1:n2, 1:n3, OV);
sum(msk_LV(:))/itpFactor*(prod(dcmInfo{1}.PixelSpacing)*dcmInfo{1}.SliceThickness)/1000

%% visualize
% LVshowresult=[LG_scale(im3d),double(im3d_bw),double(msk_LV)];
LVshowresult=[LG_scale(im3d),double(msk_LV)];
% debug view:   [permute(LG_scale(im3d),[3 2 1]),permute(double(msk_LV),[3 2 1])]



%% junk
% list_cineInds=[1:9,11,13:16]; % which of the folders listed by dir('sax*') to load
% list_cines=dir([path_case,filesep,'sax*']);
% list_cines=list_cines(list_cineInds);
