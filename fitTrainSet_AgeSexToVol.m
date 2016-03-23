%% intro

% Fits a coarse model to: age, sex, heart axis length
% 
% Expects user having dragged an existing .mat to workspace, containing:
% all_dcmInfo, all_cinedirs, all_subjdirs


%% load DICOM info .mat
nsubj=length(all_subjdirs);

all_axLen_raw=nan(nsubj,1);
all_ageYr_raw=nan(nsubj,1);
all_isFem_raw=nan(nsubj,1);

for isubj=1:nsubj
    fprintf('Subj %s:  ',all_subjdirs(isubj).name);
    % exclude: 123, 234, 277, 362, 499
    if     length(all_subjdirs(isubj).name)==1, matchStr=['00',all_subjdirs(isubj).name];
    elseif length(all_subjdirs(isubj).name)==2, matchStr=['0', all_subjdirs(isubj).name];
    else matchStr=all_subjdirs(isubj).name;
    end
    if ~isempty(strfind('123, 234, 277, 362, 499',matchStr))
        fprintf('Excluded case. Skipping. \n');
    else
        
        % pull data: axis length
        cc_slcLoc=cellfun(@(x)x.SliceLocation, all_dcmInfo{isubj});
        all_axLen_raw(isubj)=abs(max(cc_slcLoc)-min(cc_slcLoc));
        
        % pull data: age in years
        all_ageYr_raw(isubj)=getAgeYr(all_dcmInfo{isubj}{1});
        
        % pull data: gender
        all_isFem_raw(isubj)=double(~strcmpi(all_dcmInfo{isubj}{1}.PatientSex,'M'));
        fprintf('Done. \n');
    end
end

% sort data from dir() order to ascending order (as in .csv):
[~,inds_dir2asc]=sort(cellfun(@(x)str2double(x),{all_subjdirs.name}));
% test: text in cells should ascend: 1 2 3 ...
%     tmpcell={all_subjdirs(inds_dir2asc).name};
all_axLen_raw=all_axLen_raw(inds_dir2asc);
all_ageYr_raw=all_ageYr_raw(inds_dir2asc);
all_isFem_raw=all_isFem_raw(inds_dir2asc);

% remove NaN (must come after sorting)
mskNan=isnan(all_axLen_raw);
all_axLen=all_axLen_raw(~mskNan);
all_ageYr=all_ageYr_raw(~mskNan);
all_isFem=all_isFem_raw(~mskNan);
nsubj_valid=length(all_axLen);


%% load volume data
vol_raw=dlmread('/Users/liheng/Datasets/DSB16/train.csv',',',1,0); % skip row 1
vol_sys=vol_raw(~mskNan,2);
vol_dia=vol_raw(~mskNan,3);

%% test plot: axis lengthh as predictor
[all_axLen_sorted,sortInds]=sort(all_axLen);
figure;
    plot(all_axLen_sorted,vol_sys(sortInds),'ro'); hold on;
    plot(all_axLen_sorted,vol_dia(sortInds),'bo');
    title('Training dataset LV volume as predicted by axis length');
    legend('sys','dia');
    xlabel('SAX stack axis length (mm)');
    ylabel('LV volume (mL)');

%% test plot: age as predictor
[all_ageYr_sorted,sortInds]=sort(all_ageYr);
figure;
    plot(all_ageYr_sorted,vol_sys(sortInds),'ro'); hold on;
    plot(all_ageYr_sorted,vol_dia(sortInds),'bo');
    title('Training dataset LV volume as predicted by age');
    legend('sys','dia');
    xlabel('Age (yr)');
    ylabel('LV volume (mL)');

%% fit linear model
tbl = table(vol_sys,vol_dia,all_axLen,all_ageYr);
% tbl = table(vol_dia,all_axLen,all_ageYr);
tbl.isFem = nominal(all_isFem);

% fit_sys=fitlm(tbl,'vol_sys~all_axLen^2+all_ageYr^2+isFem');
% fit_dia=fitlm(tbl,'vol_dia~all_axLen^2+all_ageYr^2+isFem');
% titleStr='Training dataset: vol ~ axLen^2 + ageYr^2 + isFem';
fit_sys=fitlm(tbl,'vol_sys~all_axLen+all_ageYr+isFem');
fit_dia=fitlm(tbl,'vol_dia~all_axLen+all_ageYr+isFem');
titleStr='Training dataset: vol ~ axLen + ageYr + isFem';

%% run prediction for comparison
tbl2=table(all_axLen,all_ageYr);
tbl2.isFem=nominal(all_isFem);
vol2_sys=predict(fit_sys,tbl2);
vol2_dia=predict(fit_dia,tbl2);

%% plot prediction
figure('position',[30 30 1200 500]);
    plot(vol_sys, 'ro'); hold on;
    plot(vol2_sys,'r-'); hold on;
    plot(vol_dia, 'bo'); hold on;
    plot(vol2_dia,'b-'); hold on;
    legend('sys data','sys predicted','dia data','dia predicted');
    title(titleStr,'interpreter','none');
    xlabel('Training set case numbers');
    ylabel('LV volume (mL)');

%% run prediction on Validation or Test set
% HOW TO:
% Drag Validation set DICOM info .mat to workspace (this will overwrite
% Training set's .mat). Run "load DICOM info" cell. Then run this cell.
if false
    %% run by hand:
    tbl2=table(all_axLen,all_ageYr);
    tbl2.isFem=nominal(all_isFem);
    vol2_sys=predict(fit_sys,tbl2);
    vol2_dia=predict(fit_dia,tbl2);
    % export:
    vol2_export=[vol2_dia,vol2_sys]';
    vol2_export=vol2_export(:);
    CDFs=softstep(repmat(0:599,[nsubj*2 1]),vol2_export,4);  % debug plot:   CDFs'    
    writeCSV(CDFs,501:700);
end


%% junk

% OLD age interpretation
%         ageStr=all_dcmInfo{isubj}{1}.PatientAge;
%         all_ageYr_raw(isubj)=str2double(ageStr(1:3));
%         if strcmpi(ageStr(end),'M')
%             all_ageYr_raw(isubj)=all_ageYr_raw(isubj)/12;
%         elseif strcmpi(ageStr(end),'W')
%             all_ageYr_raw(isubj)=all_ageYr_raw(isubj)/52;
%         end


% OLD load from dicoms
% 
% all_subjdirs=dir(path_group);
% maskExclude=cellfun(@(x)strcmp(x(1),'.'),{all_subjdirs.name}) | cellfun(@(x)~x,{all_subjdirs.isdir}); % remove anything whose name starts with . or isn't a directory
% all_subjdirs=all_subjdirs(~maskExclude);
% nsubj=length(all_subjdirs);
% 
% all_axlen=nan(nsubj,1);
% all_ageyr=nan(nsubj,1);
% all_sex=nan(nsubj,1);
% 
% for isubj=1:nsubj
%     fprintf('Subj %s:  ',all_subjdirs(isubj).name);
%     % exclude: 123, 234, 277, 362, 499
%     if ~isempty(strfind('123, 234, 277, 362, 499',all_subjdirs(isubj).name))
%         fprintf('Excluded case. Skipping. \n');
%     else        
%         path_subj=[path_group,filesep,all_subjdirs(isubj).name,filesep,'study',filesep];
%         list_cine=dir([path_subj,'sax*']); % <<<< SAX filter
%         ncine=length(list_cine);
%         cc_slcLoc=nan(ncine,1);
%         for icine=1:ncine % sax_1, sax_2, ...
%             fprintf('%s, ',list_cine(icine).name);
%             path_cine=[path_subj,list_cine(icine).name,filesep];
%             [~,dcmInfo]=loadCine(path_cine,'frames',1,'print',false,'dcmInfoOnly',true);
%             
%         end        
%         fprintf('Done. \n');
%     end
% end
% 