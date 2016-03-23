
%% training set

% path_group='/Users/liheng/Datasets/DSB16/train';
% subjnums=1:500;
% 
% % Subjects to skip:
% subjnums(173)=[];
% 
% skipcine={
%     1,{'sax_57','sax_38'};
%     123,'sax_20';
% };


%% test set

path_group='/Users/liheng/Datasets/DSB16/test';
subjnums=701:1140;

skipcine={
    708,{'sax_16','sax_17','sax_18','sax_19'};
    731,'sax_10';
    746,{'sax_6','sax_7','sax_8'};
    755,{'sax_6','sax_7','sax_8','sax_9','sax_10','sax_11'};
    760,'sax_16';
    767,{'sax_13','sax_14'};
    819,'sax_17';
    831,{'sax_6','sax_7','sax_8','sax_9','sax_10','sax_11'};
    847,{'sax_9','sax_10'};
    870,{'sax_22','sax_25','sax_26','sax_28','sax_30','sax_31','sax_32','sax_33','sax_34','sax_35','sax_36','sax_37','sax_38','sax_39','sax_40'};
    1082,'sax_20';
    856,{'sax_23','sax_24','sax_25'};
    977,'sax_49';
    1079,{'sax_41','sax_40','sax_22'};
};


%% common
op.midSlcLocPct=[35 65];%[35 65]; % defines which slices are used for initial LV center finding
op.motionRegTrimSaveFrxn=0.6;
op.ClusMskMedFiltWidth=3;
op.ClusMskTrimSaveFrxn=0.9;
op.RndMskTrimThreshMax=0.95;
op.hrtMskDilateR=5;

op.crop=[0.6 0.6]; % [d1 d2], central fraction to keep

op.intensTolBld=[0.7 5]; % lower and higher bounds, in fraction of the absolute difference between blood and myo intensities.
op.landingBoxHW=15; % num of px
op.areaRatioLim_bas=2;
op.areaRatioLim_api=1.2;

op.regMskCloseRad=2;

