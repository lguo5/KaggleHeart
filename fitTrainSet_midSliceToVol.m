%% intro
% Fits a model: Kaggle training set LV sys/dias volumes from mid-slice area

%% targets

path_areas='/Users/liheng/Dropbox/KaggleBowl/midslices/areas.mat';
path_dicomInfo='/Users/liheng/Datasets/DSB16/LG/dcmInfo check - Train/set=Train N=500 frame=1.mat';

%% load
load(path_areas); % -> areas: [subj, sys, dias]
load(path_dicomInfo); % -> all_subjdirs, all_dcmInfo
vol_raw=dlmread('/Users/liheng/Datasets/DSB16/train.csv',',',1,0); % skip row 1

% path_mats='/Users/liheng/Dropbox/KaggleBowl/midslices';
% list_mats=dir([path_mats,filesep,'Patient*.mat']);
% n_mats=length(list_mats);

nsubj=size(areas,1);

%% start

% From DICOM info:
% midslc_cm2=nan(nsubj,2); % col1: sys, col2: dias.
pxmm2=nan(nsubj,1);
isFem=false(nsubj,1);
ageYr=nan(nsubj,1);
axLen=nan(nsubj,1);
% From CSV:
vol_sys=nan(nsubj,1);
vol_dia=nan(nsubj,1);

for isubj=1:nsubj
    
    % From DICOM info: 
    i_dcm=find(strcmp({all_subjdirs.name},num2str(areas(isubj,1))));
    dcmInfo=all_dcmInfo{i_dcm}{1};
    pxmm2(isubj)=prod(dcmInfo.PixelSpacing);
    isFem(isubj)=~strcmpi(dcmInfo.PatientSex,'M');
    ageYr(isubj)=getAgeYr(dcmInfo);
    
    cc_slcLoc=cellfun(@(x)x.SliceLocation, all_dcmInfo{isubj});
    axLen(isubj)=abs(max(cc_slcLoc)-min(cc_slcLoc));
    
    % From truth CSV:
    i_csv=find(vol_raw(:,1)==areas(isubj,1));
    vol_sys(isubj)=vol_raw(i_csv,2);
    vol_dia(isubj)=vol_raw(i_csv,3);
end

midslc_sys=areas(:,2).*pxmm2./100;
midslc_dia=areas(:,3).*pxmm2./100;

% hypothetical cylinder volume
cylvol_sys=midslc_sys.*axLen./10; % ml
cylvol_dia=midslc_dia.*axLen./10; % ml

%% plot
% figure;
% Vol vs. cylinder vol
%     plot(cylvol_sys,vol_sys,'ro');hold on;
%     plot(cylvol_dia,vol_dia,'bo');hold on;
%     xlabel('cylinder vol (cm^3)'); ylabel('LV Vol (ml)');
%     legend('sys','dia');
% Vol vs. cylinder vol, and age:
%     plot3(cylvol_sys,ageYr,vol_sys,'ro'); hold on;
%     plot3(cylvol_dia,ageYr,vol_dia,'bo'); hold on;
%     xlabel('cylinder vol (cm^3)'); ylabel('age (yr)'); zlabel('LV Vol (ml)');
%     legend('sys','dia');
% Vol vs. midslc, gender
%     plot3(midslc_sys,double(isFem),vol_sys,'ro'); hold on;
%     plot3(midslc_dia,double(isFem),vol_dia,'bo'); 
%     xlabel('mid-slice area (cm^2)'); ylabel('is female?'); zlabel('Vol (ml)');
%     legend('sys','dia','Location','best');
% Vol vs. midslc, age
%     plot3(midslc_sys,ageYr,vol_sys,'ro'); hold on;
%     plot3(midslc_dia,ageYr,vol_dia,'bo'); 
%     xlabel('mid-slice area (cm^2)'); ylabel('age (yr)'); zlabel('Vol (ml)');
%     legend('sys','dia','Location','best');
% Vol vs. midslc, axis length
%     plot3(midslc_sys,axLen,vol_sys,'ro'); hold on;
%     plot3(midslc_dia,axLen,vol_dia,'bo'); 
%     xlabel('mid-slice area (cm^2)'); ylabel('axis length (mm)'); zlabel('Vol (ml)');
%     legend('sys','dia','Location','best');
% Vol vs midslc, male only:
%     plot(midslc_sys(~isFem),vol_sys(~isFem),'ro'); hold on;
%     plot(midslc_dia(~isFem),vol_dia(~isFem),'bo'); 
% Vol vs midslc, Both gender:
%     plot(midslc_sys,vol_sys,'ro'); hold on;
%     plot(midslc_dia,vol_dia,'bo'); 

%% selective fit
% mskIn=ageYr<20;
% vol_sys=vol_sys(mskIn);
% vol_dia=vol_dia(mskIn);
% midslc_sys=midslc_sys(mskIn);
% midslc_dia=midslc_dia(mskIn);
% ageYr=ageYr(mskIn);


%% fit linear model
tbl = table(vol_sys,vol_dia,midslc_sys,midslc_dia,ageYr);

% titleStr='Training dataset: vol ~ midslc';
% fit_sys=fitlm(tbl,'vol_sys~midslc_sys')
% fit_dia=fitlm(tbl,'vol_dia~midslc_dia')

% titleStr='Training dataset: vol ~ midslc + ageYr';
% fit_sys=fitlm(tbl,'vol_sys~midslc_sys+ageYr')
% fit_dia=fitlm(tbl,'vol_dia~midslc_dia+ageYr')

% titleStr='Training dataset: vol ~ midslc + ageYr + midslc*ageYr';
% fit_sys=fitlm(tbl,'vol_sys~midslc_sys*ageYr')
% fit_dia=fitlm(tbl,'vol_dia~midslc_dia*ageYr')

tbl.isYouth=nominal(ageYr<22);
titleStr='Training dataset: vol ~ midslc + isYouth + midslc*isYouth';
fit_sys=fitlm(tbl,'vol_sys~midslc_sys*isYouth','CategoricalVar','isYouth')
fit_dia=fitlm(tbl,'vol_dia~midslc_dia*isYouth','CategoricalVar','isYouth')


% run prediction for comparison
tbl2_sys=table(midslc_sys,ageYr);
tbl2_dia=table(midslc_dia,ageYr);
tbl2_sys.isYouth=nominal(ageYr<20);
tbl2_dia.isYouth=nominal(ageYr<20);
vol2_sys=predict(fit_sys,tbl2_sys);
vol2_dia=predict(fit_dia,tbl2_dia);

figure;
    % data points:
    plot3(midslc_sys,ageYr,vol_sys,'ro'); hold on;
    plot3(midslc_dia,ageYr,vol_dia,'bo'); 
    % error distance as line connectors:
    plot3([midslc_sys,midslc_sys]',[ageYr,ageYr]',[vol_sys,vol2_sys]','r-'); hold on;
    plot3([midslc_dia,midslc_dia]',[ageYr,ageYr]',[vol_dia,vol2_dia]','b-'); hold on;
    % labels:
    title(titleStr,'interpreter','none');
    xlabel('mid-slice area (cm^2)'); ylabel('age (yr)'); zlabel('Vol (ml)');
    legend('sys data','dia data','sys model error','dia model error','Location','best');








