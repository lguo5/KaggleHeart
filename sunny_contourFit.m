%% intro



%% target
path_i={'/Users/liheng/Dropbox/KaggleBowl/Sunnybrook contours/training set';
        '/Users/liheng/Dropbox/KaggleBowl/Sunnybrook contours/validation set'}; % "i" for endocardial contours



%% start
ngroup=length(path_i);

all_vol_sys=cell(ngroup,1);
all_vol_dia=cell(ngroup,1);
all_midslc_mm2_sys=cell(ngroup,1);
all_midslc_mm2_dia=cell(ngroup,1);
all_slcZdif=cell(ngroup,1);

for igroup=1:ngroup
    list_cases=dir([path_i{igroup},filesep,'SC*']);
    ncase=length(list_cases);
    fprintf('Loading "%s", case:\n',pathLastPart(path_i{igroup}));
    
    % cg = current group
    cg_vol_sys=nan(ncase,1);
    cg_vol_dia=nan(ncase,1);
    cg_midslc_pos=nan(ncase,1); % interp input
    cg_midslc_mm2_sys=nan(ncase,1); % interp output
    cg_midslc_mm2_dia=nan(ncase,1); % interp output
    cg_slcZdif=nan(ncase,1);
    
    for icase=1:ncase
        %% 
        fprintf('\t%s\n',list_cases(icase).name);
        path_case=[path_i{igroup},filesep,list_cases(icase).name,filesep];
        dir_mat=dir([path_case,'*.mat']);
        load([path_case,dir_mat.name]);
        
        % cc = current case:
        cc_slcLocs=cellfun(@(x)x.SliceLocation, rec_dicomInfo);
        cc_slcThck=cellfun(@(x)x.SliceThickness,rec_dicomInfo); % mm^2
        cc_pxMM2=cellfun(@(x)prod(x.PixelSpacing),rec_dicomInfo); % mm^2
   
        cc_contrMM2=rec_contourMskNpt.*cc_pxMM2;
        cc_slcLocs_u=unique(cc_slcLocs);
        cc_slcNUnique=length(cc_slcLocs_u);
        cc_sys_mm2=nan(cc_slcNUnique,1);
        cc_dia_mm2=nan(cc_slcNUnique,1);
        cc_mskValid=false(cc_slcNUnique,1);
        % sort through each unique slice location:
        for i_slcUnique=1:cc_slcNUnique
            currAreas=cc_contrMM2(cc_slcLocs_u(i_slcUnique)==cc_slcLocs);
            if length(currAreas)==2 % need 2 values at any slice location: smaller one: sys, larger: dias
                cc_sys_mm2(i_slcUnique)=min(currAreas);
                cc_dia_mm2(i_slcUnique)=max(currAreas);
                cc_mskValid(i_slcUnique)=true;
            end
        end
        % delete unfilled spots
        cc_sys_mm2=cc_sys_mm2(cc_mskValid); 
        cc_dia_mm2=cc_dia_mm2(cc_mskValid);
        cc_slcLocs_u=cc_slcLocs_u(cc_mskValid);
        cc_slcThck=cc_slcThck(cc_mskValid);
        
        % find volume
        cg_vol_sys(icase)=sum(cc_sys_mm2.*cc_slcThck);
        cg_vol_dia(icase)=sum(cc_dia_mm2.*cc_slcThck);
        
        % find mid area
        cg_midslc_pos(icase)=0.5*(cc_slcLocs_u(1)+cc_slcLocs_u(end));
        cg_midslc_mm2_sys(icase)=interp1(cc_slcLocs_u,cc_sys_mm2,cg_midslc_pos(icase));
        cg_midslc_mm2_dia(icase)=interp1(cc_slcLocs_u,cc_dia_mm2,cg_midslc_pos(icase));
        
        % axis length
        cg_slcZdif(icase)=abs(max(cc_slcLocs)-min(cc_slcLocs));
        
%         %% plot
%         if icase==1 && igroup==1
%             hFig=figure('position',[50 50 1200 500]);
%         end
%         xplot=abs(cc_slcLocs_u); xplot=xplot-min(xplot);
%         plot(xplot, cc_sys_mm2, 's-'); hold on;
%         plot(xplot, cc_dia_mm2, 's-'); hold on;
%         % no good: plot(xplot, LG_scale(cc_sys_mm2), 's-'); hold on;
%         % plot([xplot',xplot'], [cc_sys_mm2',cc_dia_mm2'], 's-'); hold on;
        
    end % loop through cases
    
    all_vol_sys{igroup}=cg_vol_sys/1E3;
    all_vol_dia{igroup}=cg_vol_dia/1E3;
    all_midslc_mm2_sys{igroup}=cg_midslc_mm2_sys;
    all_midslc_mm2_dia{igroup}=cg_midslc_mm2_dia;
    all_slcZdif{igroup}=cg_slcZdif;
    
end % loop through groups

%% convert to linear vector form for plotting and fitting
vol_sys=cell2mat(all_vol_sys(:));
vol_dia=cell2mat(all_vol_dia(:));
midslc_cm2_sys=cell2mat(all_midslc_mm2_sys(:))./100;
midslc_cm2_dia=cell2mat(all_midslc_mm2_dia(:))./100;
slcZdif=cell2mat(all_slcZdif(:));

%% plot: vol vs. mid slice area
figure;
    plot(midslc_cm2_sys, vol_sys,'rs'); hold on;
    plot(midslc_cm2_dia, vol_dia,'bs');
    legend('sys','dia');
    title('Sunnybrook Train+Val datasets: LV vol vs. mid-slice area');
    xlabel('Mid-slice area (cm^2)');
    ylabel('LV (mL)');

%% linear model fit
tbl = table(vol_sys,vol_dia,midslc_cm2_sys,midslc_cm2_dia,slcZdif);
% fit_sys=fitlm(tbl,'vol_sys~midslc_cm2_sys^2+slcZdif^2')
% fit_dia=fitlm(tbl,'vol_dia~midslc_cm2_dia^2+slcZdif^2')
% titleStr='Sunnybrook: vol~midslc_cm2^2+slcZdif^2';

fit_sys=fitlm(tbl,'vol_sys~midslc_cm2_sys+slcZdif')
fit_dia=fitlm(tbl,'vol_dia~midslc_cm2_dia+slcZdif')
titleStr='Sunnybrook: vol~midslc_cm2+slcZdif';

% self-test model:
tbl2_sys=table(midslc_cm2_sys,slcZdif);
tbl2_dia=table(midslc_cm2_dia,slcZdif);
vol2_sys=predict(fit_sys,tbl2_sys);
vol2_dia=predict(fit_dia,tbl2_dia);

% plot
figure('position',[30 30 1200 500]);
    plot(vol_sys, 'ro'); hold on;
    plot(vol2_sys,'r-'); hold on;
    plot(vol_dia, 'bo'); hold on;
    plot(vol2_dia,'b-'); hold on;
    legend('sys data','sys predicted','dia data','dia predicted');
    title(titleStr,'interpreter','none');
    xlabel('Nth case');
    ylabel('LV volume (mL)');


%% JUNK old fit
% % dependent on mid-slice area, axis length:
% X_sys=[cell2mat(all_midslc_mm2_sys(:))./100,cell2mat(all_slcZdif(:))];
% X_dia=[cell2mat(all_midslc_mm2_dia(:))./100,cell2mat(all_slcZdif(:))];
% T=[0 0 0; 1 0 0; 0 1 0;];
% 
% % dependent on mid-slice area
% X_sys=cell2mat(all_midslc_mm2_sys(:));
% X_dia=cell2mat(all_midslc_mm2_dia(:));
% T=[0 0; 1 0];
% 
% fit_sys=fitlm(X_sys,cell2mat(all_vol_sys(:)),T)
% fit_dia=fitlm(X_dia,cell2mat(all_vol_dia(:)),T)
% 
% dependent on mid-slice area, axis length:
% 
%                    Estimate       SE        tStat       pValue  
%                    ________    _________    ______    __________
% 
%     (Intercept)     -65.613       19.141    -3.428     0.0019649
%     x1             0.065669    0.0034463    19.055    3.4617e-17
%     x2               1.1983      0.38803    3.0881     0.0046247
% Number of observations: 30, Error degrees of freedom: 27
% Root Mean Squared Error: 17
% R-squared: 0.949,  Adjusted R-Squared 0.945
% F-statistic vs. constant model: 250, p-value = 3.75e-18
% 
%                    Estimate       SE         tStat       pValue  
%                    ________    _________    _______    __________
% 
%     (Intercept)      -107.4       23.151    -4.6393    8.0282e-05
%     x1             0.077012    0.0055949     13.765    1.0107e-13
%     x2               1.4598      0.49292     2.9615     0.0063132
% Number of observations: 30, Error degrees of freedom: 27
% Root Mean Squared Error: 20.6
% R-squared: 0.92,  Adjusted R-Squared 0.914
% F-statistic vs. constant model: 154, p-value = 1.66e-15