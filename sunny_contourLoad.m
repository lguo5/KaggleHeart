%% intro
% 
% BEHAVIOR:
%     For each case in contour data (_con):
%         Finds corresponding folder in DICOM data (_dcm)
%         Makes a folder in current folder (w caseName)
%         For each contour text file:
%             Finds matching DICOM in corresponding folder
%             Show DICOM w contour overlaid, saves to the newly created case folder
%             Records DICOM info
%         Saves contour info and DICOM info to the newly created case folder
% 
% 
saveJPG_DICOMwContour=true;


%% target selection

% % Sunnybrook Cardiac MR Database ContoursPart3
% path_group_con='/Users/liheng/Datasets/Sunnybrook/Sunnybrook Cardiac MR Database ContoursPart3/TrainingDataContours';
% path_group_dcm='/Users/liheng/Datasets/Sunnybrook/challenge_training';
% filt_con='i'; % filter to only allow endo contour: as in IM-0001-0068-iXXXX.txt
% subjtags={ % name of folders in dicom data corresponding to contour folders (but without padded 0s), and series number of scans that are actually contoured.
%     'SC-HF-I-1', '0004';
%     'SC-HF-I-2', '0106';
%     'SC-HF-I-4', '0116';
%     'SC-HF-I-40', '0134';
%     'SC-HF-NI-3', '0379';
%     'SC-HF-NI-4', '0501';
%     'SC-HF-NI-34', '0446';
%     'SC-HF-NI-36', '0474';
%     'SC-HYP-1', '0550';
%     'SC-HYP-3', '0650';
%     'SC-HYP-38', '0734';
%     'SC-HYP-40', '0755';
%     'SC-N-2', '0898';
%     'SC-N-3', '0915';
%     'SC-N-40', '0944'};

% % Sunnybrook Cardiac MR Database ContoursPart2
% path_group_con='/Users/liheng/Datasets/Sunnybrook/Sunnybrook Cardiac MR Database ContoursPart2/ValidationDataContours';
% path_group_dcm='/Users/liheng/Datasets/Sunnybrook/challenge_validation';
% filt_con='i'; % filter to only allow endo contour: as in IM-0001-0068-iXXXX.txt
% subjtags={ % name of folders in dicom data corresponding to contour folders (but without padded 0s), and series number of scans that are actually contoured.
%     'SC-HF-I-5','0156'; % 0174
%     'SC-HF-I-6','0180'; 
%     'SC-HF-I-7','0209'; % 217
%     'SC-HF-I-8','0226'; % 
%     'SC-HF-NI-7','0523'
%     'SC-HF-NI-11','0270'
%     'SC-HF-NI-31','0401'
%     'SC-HF-NI-33','0424' % 428
%     'SC-HYP-6','0767' % 768 770 771
%     'SC-HYP-7','0007' % 0011 
%     'SC-HYP-8','0796' % 797 
%     'SC-HYP-37','0702' % 714
%     'SC-N-5','0963' % 965
%     'SC-N-6','0984' % 985
%     'SC-N-7','1009' % 1013
% };


%% start
list_subjs_con=dir([path_group_con,filesep,'SC-*']);
list_subjs_dcm=dir([path_group_dcm,filesep,'SC-*']);
n_subj_con=length(list_subjs_con);

for i_subj_con=1:n_subj_con
    %% run by hand
    caseName_con =list_subjs_con(i_subj_con).name;
    caseName_dcm=stripNonEnding0(caseName_con);
    caseMask_dcm=strcmp(subjtags(:,1),caseName_dcm);
    if sum(caseMask_dcm)==0
        warning('Cound not find DICOM folder for contour case "%s". Skipping this contour case...\n',caseName_con);
    elseif sum(caseMask_dcm)>1
        warning('Found more than 1 DICOM folders for contour case "%s". Skipping this contour case...\n',caseName_con);
    else
        fprintf('Found 1 DICOM folder for contour case "%s"...\n',caseName_con);
        mkdir(caseName_con); % this does not overwrite existing folder
        filt_dcm=subjtags{caseMask_dcm,2}; % e.g. '0004'
        
        %% looping through contour files
        path_case_con=[path_group_con,filesep,caseName_con,filesep,'contours-manual/IRCCI-expert'];
        list_cons=dir(path_case_con);
        msk=cellfun(@(x)~isempty(regexp(x,['IM-0001-\d\d\d\d-',filt_con])),{list_cons.name});
        list_cons=list_cons(msk);
        n_con=length(list_cons);
        
        rec_init=true; %
        for i_con=1:n_con
            filt_frm=list_cons(i_con).name(9:12); % e.g. '0048'
            path_1dcm=[path_group_dcm,filesep,caseName_dcm,filesep,'IM-',filt_dcm,'-',filt_frm,'.dcm'];
            list_1dcm=dir(path_1dcm);
            if isempty(list_1dcm)
                warning('\tDid not find matching DICOM for:  "%s"... Skipping contour file.\n',list_cons(i_con).name);
            else
                matchDicom=dicomread(path_1dcm);
                matchDicomInfo=dicominfo(path_1dcm);
                % load contour and convert to mask:
                contourVtx=dlmread([path_case_con,filesep,list_cons(i_con).name]);
                contourMsk=poly2mask(contourVtx(:,1),contourVtx(:,2),double(matchDicomInfo.Rows),double(matchDicomInfo.Columns));
                %% show/save image + contour
                if saveJPG_DICOMwContour
                    hFig=figure('paperpositionmode','auto');
                        imshow(LG_normalize(double(matchDicom)),'Border','tight'); hold on;
                        plot(contourVtx(:,1),contourVtx(:,2));
                        print('-djpeg',sprintf('%s%s%s, %s',caseName_con,filesep,list_1dcm.name(1:end-4),list_cons(i_con).name(1:end-4)),'-r300');
                        close(hFig);
                end
                if rec_init % use this flag instead of ind==1, lest no DICOM found for ind==1
                    rec_contourVtx=cell(n_con,1);
                    rec_contourMsk=cell(n_con,1);
                    rec_contourMskNpt=nan(n_con,1);
                    rec_contourFileName=cell(n_con,1);
                    rec_dicomFileName=cell(n_con,1);
                    rec_dicomInfo=cell(n_con,1);
                end
                rec_init=false;
                rec_contourVtx{i_con}=contourVtx;
                rec_contourMsk{i_con}=contourMsk;
                rec_contourMskNpt(i_con)=sum(contourMsk(:));
                rec_contourFileName{i_con}=list_cons(i_con).name;
                rec_dicomFileName{i_con}=list_1dcm.name;
                rec_dicomInfo{i_con}=matchDicomInfo;
            end % found matching DICOM or not
        end % contour file loop
        save(sprintf('%s%s%s contours',caseName_con,filesep,caseName_con),'rec*');
        
        %% make/save graph to folder
        % You can drag a saved .mat to workspace and run this cell
        msk=cellfun(@(x)~isempty(x),rec_dicomInfo); % to skip any contour w/o a matching DICOM
        slcLocs=cellfun(@(x)x.SliceLocation,rec_dicomInfo(msk)); % cell2mat({rec_dicomInfo.SliceLocation});
        pxAreas=cellfun(@(x)prod(x.PixelSpacing),rec_dicomInfo(msk));
        contourAreas=rec_contourMskNpt(msk).*pxAreas;

        hFig=figure('position', [22 50 560 240], 'paperpositionmode', 'auto');
            plot(slcLocs,contourAreas,'s');
            title(caseName_con,'interpreter','none');
            xlabel('Slice Location (mm)');
            ylabel('Contour Area (mm^2)');
            axis tight;
            print('-dpng',sprintf('%s%s%s contour area',caseName_con,filesep,caseName_con),'-r80');
            close(hFig);

        fprintf('Done. Found matching DICOMs for %d of %d contours.\n',sum(msk),length(msk));
        
    end % found matching DICOM folder or not
end % subject loop

%%
% matchTb=table({'SC-HF-I-1';	'SC-HF-I-2';	'SC-HF-I-4';	'SC-HF-I-40';	'SC-HF-NI-3';	'SC-HF-NI-4';	'SC-HF-NI-34';	'SC-HF-NI-36';	'SC-HYP-1';	'SC-HYP-3';	'SC-HYP-38';	'SC-HYP-40';	'SC-N-2';	'SC-N-3';	'SC-N-40'},...
%     {'0004';	'0106';	'0116';	'0134';	'0379';	'0501';	'0446';	'0474';	'0550';	'0650';	'0734';	'0755';	'0898';	'0915';	'0944'},...
%     'VariableNames',{'subj','tag'});
% 
% matchTb=table({'0004';	'0106';	'0116';	'0134';	'0379';	'0501';	'0446';	'0474';	'0550';	'0650';	'0734';	'0755';	'0898';	'0915';	'0944'},...
%     'RowNames',{'SC-HF-I-1';	'SC-HF-I-2';	'SC-HF-I-4';	'SC-HF-I-40';	'SC-HF-NI-3';	'SC-HF-NI-4';	'SC-HF-NI-34';	'SC-HF-NI-36';	'SC-HYP-1';	'SC-HYP-3';	'SC-HYP-38';	'SC-HYP-40';	'SC-N-2';	'SC-N-3';	'SC-N-40'},...
%     'VariableNames',{'tag'});

