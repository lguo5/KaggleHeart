
%% target selection
runGroup_ops();

%% run control
saveFolder_top=['regGroTest_',pathLastPart(path_group)];
% saveFolder_top=['midSlcGrpTest_',pathLastPart(path_group)];
savename_log='log.txt';
savename_vol='vol.csv';

%% start
mkdir(saveFolder_top);
fid=fopen([saveFolder_top,filesep,savename_log],'a'); % <<<<< append?
fid2=fopen([saveFolder_top,filesep,savename_vol],'a'); % <<<<< append?
% do not write header row due to append:
%     fprintf(fid2,'Id,Systole,Diastole\n');

nsubj=length(subjnums);
% 123: all-in-1 case
% 337: PE size change  
% 442: Last 5 cines transposed
for isubj=832:1140 % subjnums  % isubj is directly translated to sax_ name
    
    logstr=sprintf('Subj %d:\n',isubj);    fprintf(logstr);
    if exist('fid','var'), fprintf(fid,logstr); end
    
    % List all SAX folders
    name_subj=num2str(isubj);
    path_study=[path_group,filesep,name_subj,filesep,'study',filesep];
    list_cines=dir([path_study,'sax*']);
    
    % Remove those hand-marked "to skip"
%     skipname=skipcine_cinename(skipcine_subjname==isubj);
    skipname=skipcine(cell2mat(skipcine(:,1))==isubj,2); % must use () here to extract to 1x1 cell, then extract w {1} later
    if ~isempty(skipname)
        skiplist=cellfun(@(x)any(strcmpi(skipname{1},x)),{list_cines.name});
        list_cines(skiplist)=[];  % test >>>>    {list_cines.name}     
    end
    ncine=length(list_cines);
    
    % Count num of DICOMs in each SAX folder:
    list_ndcm=zeros(ncine,1); % listed in dir() order
    for icine=1:ncine
        list_dcm=dir([path_study,list_cines(icine).name,filesep,'*.dcm']);
        list_ndcm(icine)=length(list_dcm);
    end
    
    %% load cines, sort wrt slice locations
    %     __      __      __ __  __ ___ 
    % |  /  \ /\ |  \    (_ /  \|__) |  
    % |__\__//--\|__/    __)\__/| \  |  
    % Common output: 
    %    nslice: num of slices for subject
    %    slcs_im4d: 1d cell x 3d img
    %    slcs_locs: sorted neg to pos
    %    slcs_nfr: num of frames per slice
    %    slcs_pxsize: [w h] for slices listed in "slicelocs"
    %    slcs_thk: thicknesses for slices listed in "slicelocs"
    
    % All-in-one folder: only loads the FIRST all-in-one folder listed
    % NO 2 cines should have same slice location in all-in-one dicoms.
    allInOne=any(list_ndcm>200);
    if allInOne
        allInOne_folder=list_cines(find(list_ndcm>200,1,'first')).name;
        path_cine=[path_study,allInOne_folder];
        logstr=sprintf('  Loading ALL-IN-ONE "%s": ',allInOne_folder);    fprintf(logstr);
        if exist('fid','var'), fprintf(fid,logstr); end
        [load_im, load_dcmInfo]=loadCine(path_cine);
        load_slicelocs=cellfun(@(x)x.SliceLocation,load_dcmInfo); % load_dcmInfo is cell
        % load_slicelocs=cell2mat({load_dcmInfo.SliceLocation}); % load_dcmInfo is struct array
        slcs_loc=unique(load_slicelocs); slcs_loc=slcs_loc(:);
        nslice=length(slcs_loc);
        % allocate
        slcs_im4d=cell(nslice,1);
        slcs_nfr=zeros(nslice,1);
        slcs_pxsize=zeros(nslice,2);
        slcs_thk=zeros(nslice,1);
        slcs_nd2=zeros(nslice,1);
        slcs_nd1=zeros(nslice,1);
        % crop
        n1o2=size(load_im,1)/2;
        n2o2=size(load_im,2)/2;
        for islice=1:nslice
            msk=load_slicelocs==slcs_loc(islice);
            slcs_im4d{islice}=load_im((n1o2-round(n1o2*op.crop(1))+1):(n1o2+round(n1o2*op.crop(1))),...
                                      (n2o2-round(n2o2*op.crop(2))+1):(n2o2+round(n2o2*op.crop(2))),msk);
            slcs_nfr(islice)=size(slcs_im4d{islice},3);
            slcs_pxsize(islice,1)=cellfun(@(x)x.PixelSpacing(1),load_dcmInfo(find(msk,1,'first'))); % load_dcmInfo is cell
            slcs_pxsize(islice,2)=cellfun(@(x)x.PixelSpacing(2),load_dcmInfo(find(msk,1,'first'))); % load_dcmInfo is cell
            slcs_thk(islice)=cellfun(@(x)x.SliceThickness,load_dcmInfo(find(msk,1,'first')));
            % slcs_pxsize(islice,:)=load_dcmInfo(find(msk,1,'first')).PixelSpacing; % load_dcmInfo is struct array
            % slcs_thk(islice)=load_dcmInfo(find(msk,1,'first')).SliceThickness; % load_dcmInfo is struct array
            slcs_nd1(islice)=size(slcs_im4d{islice},1);
            slcs_nd2(islice)=size(slcs_im4d{islice},2);
        end
        slcs_nd1_raw=slcs_nd1; % all-in-one mode needs no transposing
        slcs_nd2_raw=slcs_nd2; % all-in-one mode needs no transposing
        
    % Slice-per-folder mode:
    else
        
        nslice=ncine;
        logstr=sprintf('  Loading by folder: ');    fprintf(logstr);
        if exist('fid','var'), fprintf(fid,logstr); end
        
        % allocate --------
        slcs_im4d=cell(nslice,1);
        slcs_loc=zeros(nslice,1);
        slcs_nfr=zeros(nslice,1);
        slcs_pxsize=zeros(nslice,2);
        slcs_thk=zeros(nslice,1);
        slcs_nd2=zeros(nslice,1);
        slcs_nd1=zeros(nslice,1);
        % -----------------
        
        for islice=1:nslice
            logstr=sprintf('%s=',list_cines(islice).name);     fprintf(logstr);
            if exist('fid','var'), fprintf(fid,logstr); end
            [im3d,dcmInfo]=loadCine([path_study,list_cines(islice).name],'print',false,'dcmInfo1stFrOnly',true);
            % crop (move to after transpose?)
            n1o2=size(im3d,1)/2;
            n2o2=size(im3d,2)/2;
            slcs_im4d{islice}=im3d((n1o2-round(n1o2*op.crop(1))+1):(n1o2+round(n1o2*op.crop(1))),...
                                   (n2o2-round(n2o2*op.crop(2))+1):(n2o2+round(n2o2*op.crop(2))),:);
            slcs_nfr(islice)=size(slcs_im4d{islice},3);
            slcs_nd1(islice)=size(slcs_im4d{islice},1);
            slcs_nd2(islice)=size(slcs_im4d{islice},2);
            slcs_loc(islice)=dcmInfo{1}.SliceLocation;
            slcs_pxsize(islice,:)=dcmInfo{1}.PixelSpacing;
            slcs_thk(islice)=dcmInfo{1}.SliceThickness;
            logstr=sprintf('%dfr  ',slcs_nfr(islice));    fprintf(logstr);
            if exist('fid','var'), fprintf(fid,logstr); end
        end
        if exist('fid','var'), fprintf(fid,'\n'); end
        fprintf('\n');
        
        % sort slice loc neg to pos; REPEATED SLICE LOC POSSIBLE and is OK
        [slcs_loc,sortInds]=sort(slcs_loc);
        slcs_im4d=slcs_im4d(sortInds);
        slcs_nfr=slcs_nfr(sortInds);
        slcs_pxsize=slcs_pxsize(sortInds,:);
        slcs_thk=slcs_thk(sortInds);
        slcs_nd1=slcs_nd1(sortInds);
        slcs_nd2=slcs_nd2(sortInds);
        list_cines=list_cines(sortInds); % test >>>>    {list_cines.name}     
        
        %% handle transpose
        % ___ __          __ __  __  __ __ 
        %  | |__) /\ |\ |(_ |__)/  \(_ |_  
        %  | | \ /--\| \|__)|   \__/__)|__ 
        slcs_d2od1=round(slcs_nd2./slcs_nd1,4);
        d2od1_mode=mode(slcs_d2od1);
        msk_slcAdjust=slcs_d2od1~=d2od1_mode;
        inds_slcAdjust=find(msk_slcAdjust);
        inds_slcNoAdjust=find(~msk_slcAdjust);
        slcs_d1od2=slcs_nd1./slcs_nd2; % the inverse
        slcs_nd1_raw=slcs_nd1; % make copy of old dims
        slcs_nd2_raw=slcs_nd2; % make copy of old dims
        if ~isempty(inds_slcAdjust)
            logstr=sprintf('  Transposing (sorted slice inds): %s: ',mat2str(inds_slcAdjust'));     fprintf(logstr);
            if exist('fid','var'), fprintf(fid,logstr); end
            for islice=inds_slcAdjust' % support is all cines
                logstr=sprintf('%s ',list_cines(islice).name);     fprintf(logstr);
                if exist('fid','var'), fprintf(fid,logstr); end
                % if inverse ratio is closer to majority ratio, rotate 90 CW:
                if abs(slcs_d1od2(islice)-d2od1_mode) < abs(slcs_d2od1(islice)-d2od1_mode) 
                    [~,isort]=min(abs(inds_slcNoAdjust-islice)); % find nearest non-adj slice as reference
                    slcs_im4d{islice}=imgRot90CWorCCW(slcs_im4d{islice},slcs_im4d{inds_slcNoAdjust(isort)});
                    % slcs_im4d{islice}=flip(permute(slcs_im4d{islice},[2 1 3]),2);
                    slcs_nd1(islice)=size(slcs_im4d{islice},1);
                    slcs_nd2(islice)=size(slcs_im4d{islice},2);
                    slcs_pxsize(islice,:)=flip(slcs_pxsize(islice,:),2);
                    slcs_nd1(islice)=slcs_nd2_raw(islice); % note LHS is 1, RHS is 2
                    slcs_nd2(islice)=slcs_nd1_raw(islice); % note LHS is 2, RHS is 1
                end
            end
            if exist('fid','var'), fprintf(fid,'\n'); end
            fprintf('\n');
        end % if transposing any image        
    end
    
%     %% print sorted folder and slice loc:
%     fprintf('  Sorted folder, sliceloc:  ');
%     for islice=1:nslice
%         fprintf('%s=%.3f  ',list_cines(islice).name, slcs_loc(islice));
%     end
%     fprintf('\n');
    
    %% decide which direction is apical
    slcs_frmNorm=cell(nslice,1);
    slcs_norm=nan(nslice,1);
    % slcs_frmNormRange=nan(nslice,1);
    % slcs_frmNormStdev=nan(nslice,1);
    slcs_frmTdifTotal=nan(nslice,1);
    for islice=1:nslice
        curr_frmNorm=nan(slcs_nfr(islice),1);
        curr_imReshape=reshape(slcs_im4d{islice},[slcs_nd2(islice)*slcs_nd1(islice),slcs_nfr(islice)]);
        for ifr=1:slcs_nfr(islice)
            curr_frmNorm(ifr)=norm(curr_imReshape(:,ifr));
        end
        slcs_frmNorm{islice}=curr_frmNorm;
        slcs_norm(islice)=norm(curr_frmNorm);
        % slcs_frmNormRange(islice)=max(curr_frmNorm)-min(curr_frmNorm);
        % slcs_frmNormStdev(islice)=std(curr_frmNorm);
        slcs_frmTdifTotal(islice)=sum(sum(abs(diff(curr_imReshape,[],2)),1));
    end
    % Find whether direction of increasing slice index is basal->apical:
    dirB2A=-sign(sum(diff(slcs_frmTdifTotal))); % if data points are majority inreasing
    
    % apical2basal=TRUE means:  sliceloc neg->pos  =  apical->basal
    % apical2basal=sum(diff(slcs_frmTdifTotal)>0)>0; % if data points are majority inreasing
    
    
    %% print sorted folders and their properties
    if allInOne, slcs_foldername=repmat({allInOne_folder},[nslice,1]);
    else         slcs_foldername={list_cines.name}';
    end
    slcs_summary=table((1:nslice)',slcs_foldername,round(slcs_loc,6),slcs_nd1_raw,slcs_nd2_raw,slcs_nd1,slcs_nd2,round(slcs_norm./slcs_nd2(islice)./slcs_nd1(islice),3),...
        'VariableNames',{'SlcInd','Folder','SlcLoc','nd1Raw','nd2Raw','nd1','nd2','pxnorm'});
    disp(slcs_summary);
    writetable(slcs_summary,sprintf('%s%s%04d slcTable.txt',saveFolder_top,filesep,isubj),'WriteRowNames',true,'Delimiter','tab');
    
    %% Find LV center map for ALL slices (but not all are used to find actual LV center)    
    %            __ __    ___ __ __              __  __ 
    % | \  /    /  |_ |\ | | |_ |__)    |\/| /\ |__)(_  
    % |__\/     \__|__| \| | |__| \     |  |/--\|   __) 
    slcs_ctrMap=cell(nslice,1);
    slcs_intens=nan(nslice,2); % col1: myocardium. col2: blood
    % slcs_ctrPos=nan(nslice,2); % dim1 fraction, dim2 fraction. Fraction is ind1/nd1, ind2/nd2, from upper left corner of image.
    
    % print
    logstr=sprintf('  Finding LV centers map: ');  fprintf(logstr);
    if exist('fid','var'), fprintf(fid,logstr); end    
    
    for islice=1:nslice%midSlcs_inds
        % print
        logstr=sprintf('%02d=%s ',islice,slcs_foldername{islice});  fprintf(logstr);
        if exist('fid','var'), fprintf(fid,logstr); end
        % LV center map:
        [slcs_ctrMap{islice}, slcs_intens(islice,1), slcs_intens(islice,2), reportImg]=...
            findLVctr_cluster(slcs_im4d{islice},...
                'ClusMskMedFiltWidth',op.ClusMskMedFiltWidth,...
                'ClusMskTrimSaveFrxn',op.ClusMskTrimSaveFrxn,...
                'RndMskTrimThreshMax',op.RndMskTrimThreshMax,...
                'motionRegTrimSaveFrxn',op.motionRegTrimSaveFrxn,...
                'hrtMskDilateR',op.hrtMskDilateR);
        % save report image
        % imwrite(reportImg,sprintf('%s%s%04d mid slc%02d(%s).jpg',saveFolder_top,filesep,isubj,islice,slcs_foldername{islice}),'jpg','quality',80);
    end
    fprintf('\n');
    if exist('fid','var'), fprintf(fid,'\n'); end
    
    %% find composite LV center from middle slice group, apply to ALL slices
    % output: midSlcs_c1c2 [1 1], midSlcs_ctrMapMeanMarked [n1min n2min]
    
    midSlcs_msk=slcs_loc>=prctile(slcs_loc,op.midSlcLocPct(1)) & slcs_loc<=prctile(slcs_loc,op.midSlcLocPct(2));
    midSlcs_inds=find(midSlcs_msk);
    midSlcs_inds=midSlcs_inds(:)'; % force to row vec
    midSlcs_n=length(midSlcs_inds);
    
    n1min=min(slcs_nd1);
    n2min=min(slcs_nd2);
    midSlcs_ctrMapMean=zeros(n1min,n2min);
    for islice=midSlcs_inds
        midSlcs_ctrMapMean=midSlcs_ctrMapMean+imgCropOrExpand(slcs_ctrMap{islice},n1min,n2min);
    end
    midSlcs_ctrMapMean=midSlcs_ctrMapMean./midSlcs_n;
    % composite center based on middle slice group
    [midSlcs_c1c2,midSlcs_ctrMapMeanMarked]=imgFindMaxRegCtr(midSlcs_ctrMapMean,0.9);
    
    % export MP4
    %      __           __          __  
    % |\/||__) |__|    |  \/  \|\/||__) 
    % |  ||       |    |__/\__/|  ||    
    reportMP4Ctr=zeros(n1min,n2min*2,nslice);
    for islice=1:nslice
        tmpImg=imgCropOrExpand(LG_scale(slcs_im4d{islice}(:,:,1)),n1min,n2min);
        tmpImg(midSlcs_c1c2(1),:)=1;
        tmpImg(:,midSlcs_c1c2(2))=1;
        reportMP4Ctr(:,:,islice)=[tmpImg,midSlcs_ctrMapMeanMarked];
    end
    % dumpMP4_img3D(reportMP4Ctr,sprintf('%s%s%04d ctrs',saveFolder_top,filesep,isubj),10);
    
    %% find SYS and DIAS frames
    %  __      __      __           __    __     __ __          __ __ 
    % |_ ||\ ||  \    |  \| /\     (_ \_/(_     |_ |__) /\ |\/||_ (_  
    % |  || \||__/    |__/|/--\    __) | __)    |  | \ /--\|  ||____) 
    % (unified for ALL cines, since basal sys frames can actually have highest intensity)
    slcs_ifr_dia=nan(nslice,1); % midSlcs_indFrmDia=nan(midSlcs_n,1);
    slcs_ifr_sys=nan(nslice,1); % midSlcs_indFrmSys=nan(midSlcs_n,1);
    for islice=midSlcs_inds
        % find DIAS frames
        [~,slcs_ifr_dia(islice)]=max(slcs_frmNorm{islice});
        % If dias frame is later in the cine (e.g. frame 28), shift up by a period so
        % it's near frame 0 (e.g. -2), so that it falls near other detected dias
        % frames (e.g. 1 or 2) for meaningful median/mean calculation:
        if slcs_ifr_dia(islice)>slcs_nfr(islice)/2,
            slcs_ifr_dia(islice)=slcs_ifr_dia(islice)-slcs_nfr(islice);
        end
        % find SYS frames:
        % No need to consider 1-period shift, since sys is usually in middle of cardiac cycle.
        [~,slcs_ifr_sys(islice)]=min(slcs_frmNorm{islice});
    end
    subj_iFr_sys=round(median(slcs_ifr_sys(midSlcs_inds)));
    subj_iFr_dia=round(median(slcs_ifr_dia(midSlcs_inds)));
    subj_iFr_dia=min(subj_iFr_dia,slcs_nfr(islice));% wrap to max nfr for the cine, for unlikely case where each slice's cine has a diff nfr.
    if subj_iFr_dia<1, subj_iFr_dia=subj_iFr_dia+slcs_nfr(islice); end
    
    % print
    logstr=sprintf('  DIA frame inds from middle slices %s: %s. Final choice: %d.\n',mat2str(midSlcs_inds),mat2str(slcs_ifr_dia(midSlcs_inds)),subj_iFr_dia);  fprintf(logstr);
    if exist('fid','var'), fprintf(fid,logstr); end
    logstr=sprintf('  SYS frame inds from middle slices %s: %s. Final choice: %d.\n',mat2str(midSlcs_inds),mat2str(slcs_ifr_sys(midSlcs_inds)),subj_iFr_sys);  fprintf(logstr);
    if exist('fid','var'), fprintf(fid,logstr); end
    
    %% find area of mid slice group
    %  __      __         __  __     
    % |_ ||\ ||  \    /\ |__)|_  /\  
    % |  || \||__/   /--\| \ |__/--\ 
    slcs_LVnpx_raw_dia=nan(nslice,1);
    slcs_LVnpx_raw_sys=nan(nslice,1);
    slcs_LVmsk_raw_dia=false(n1min,n2min,nslice);%slcs_LVmsk_dia=cell(nslice,1);
    slcs_LVmsk_raw_sys=false(n1min,n2min,nslice);%slcs_LVmsk_sys=cell(nslice,1);
    slcs_LVmskMk_dia=false(n1min,n2min,nslice);
    slcs_LVmskMk_sys=false(n1min,n2min,nslice);
    slcs_landing_dia=nan(nslice,2); %[d1 d2]
    slcs_landing_sys=nan(nslice,2); %[d1 d2]
    
    reportMP4RegGroRaw=zeros(n1min*2,n2min*2,nslice); 
    % [dia frame, dia region growth;
    %  sys frame, sys region growth]
    
    % print
    logstr=sprintf('  Finding LV areas:\n');  fprintf(logstr);
    if exist('fid','var'), fprintf(fid,logstr); end
    
    %
    %  __  __ __   __          __  __  __     ___     
    % |__)|_ / _ |/  \|\ |    / _ |__)/  \|  | | |__| 
    % | \ |__\__)|\__/| \|    \__)| \ \__/|/\| | |  |     
    for islice=1:nslice%midSlcs_inds
        %% find intensity bounds of blood
        intensBldMyoGap=abs(slcs_intens(islice,2)-slcs_intens(islice,1));
        intensBndBld=[slcs_intens(islice,2)-intensBldMyoGap.*op.intensTolBld(1),...
                      slcs_intens(islice,2)+intensBldMyoGap.*op.intensTolBld(2)];
        tmpImg3d=LG_scale(slcs_im4d{islice}); % imitates what happened in findLVctr_cluster()
        
        % print
        logstr=sprintf('    Slice%02d: BloodIntens[lo,ctr,hi]=[%.3f,%.3f,%.3f].\n',islice,intensBndBld(1),slcs_intens(islice,2),intensBndBld(2));  fprintf(logstr);
        if exist('fid','var'), fprintf(fid,logstr); end
        
        %% DIA
        tmpImg=tmpImg3d(:,:,subj_iFr_dia);
        [tmpMsk,slcs_LVnpx_raw_dia(islice),slcs_landing_dia(islice,1),slcs_landing_dia(islice,2)]=...
            findLVarea_RegGro(tmpImg,midSlcs_c1c2(1),midSlcs_c1c2(2),...
                slcs_intens(islice,2),intensBndBld,op.landingBoxHW,...
                'pxval',intensBndBld,'intensCtr',slcs_intens(islice,2));
            
        % Close little holes
        if op.regMskCloseRad>0
            tmpMsk=imclose(tmpMsk,strel('disk',op.regMskCloseRad));
        end
        
        slcs_LVmsk_raw_dia(:,:,islice)=imgCropOrExpand(tmpMsk,n1min,n2min); %slcs_LVmsk_dia{islice}
        
        % mark initial and final landing points:
        tmpMskMk=imgMarkCrossHair(tmpMsk,midSlcs_c1c2,'dashed'); % mark initial with dashed lines
        tmpMskMk=imgMarkCrossHair(tmpMskMk,slcs_landing_dia(islice,:),'solid'); % mark final with solid lines
        slcs_LVmskMk_dia(:,:,islice)=imgCropOrExpand(tmpMskMk,n1min,n2min);
        
        % print
        logstr=sprintf('      DIA landing: initial=%s final=%s\n', mat2str(midSlcs_c1c2), mat2str(slcs_landing_dia(islice,:)));  fprintf(logstr);
        if exist('fid','var'), fprintf(fid,logstr); end
        
        % export MP4
        reportMP4RegGroRaw(1:n1min,:,islice)=[imgCropOrExpand(tmpImg,n1min,n2min), slcs_LVmskMk_dia(:,:,islice)];
        
        
        %% SYS
        tmpImg=tmpImg3d(:,:,subj_iFr_sys);
        [tmpMsk,slcs_LVnpx_raw_sys(islice),slcs_landing_sys(islice,1),slcs_landing_sys(islice,2)]=...
            findLVarea_RegGro(tmpImg,midSlcs_c1c2(1),midSlcs_c1c2(2),...
                slcs_intens(islice,2),intensBndBld,op.landingBoxHW,...
                'pxval',intensBndBld,'intensCtr',slcs_intens(islice,2));
            
        % Close little holes
        if op.regMskCloseRad>0
            tmpMsk=imclose(tmpMsk,strel('disk',op.regMskCloseRad));
        end
            
        slcs_LVmsk_raw_sys(:,:,islice)=imgCropOrExpand(tmpMsk,n1min,n2min); %slcs_LVmsk_sys{islice}
        
        % mark initial and final landing points:
        tmpMskMk=imgMarkCrossHair(tmpMsk,midSlcs_c1c2,'dashed'); % mark initial with dashed lines
        tmpMskMk=imgMarkCrossHair(tmpMskMk,slcs_landing_sys(islice,:),'solid'); % mark final with solid lines
        slcs_LVmskMk_sys(:,:,islice)=imgCropOrExpand(tmpMskMk,n1min,n2min);
        
        % print
        logstr=sprintf('      SYS landing: initial=%s final=%s\n', mat2str(midSlcs_c1c2), mat2str(slcs_landing_sys(islice,:)));  fprintf(logstr);
        if exist('fid','var'), fprintf(fid,logstr); end
        
        % export MP4
        reportMP4RegGroRaw(n1min+1:end,:,islice)=[imgCropOrExpand(tmpImg,n1min,n2min), slcs_LVmskMk_sys(:,:,islice)];
        
    end
    % Dump MP4
    % dumpMP4_img3D(reportMP4RegGroRaw,sprintf('%s%s%04d regs raw',saveFolder_top,filesep,isubj),10);
    
    %% correct region growth
    %  __ __  __  __  __ _____     __  __ __   __          __  __  __     ___     
    % /  /  \|__)|__)|_ /   |     |__)|_ / _ |/  \|\ |    / _ |__)/  \|  | | |__| 
    % \__\__/| \ | \ |__\__ |     | \ |__\__)|\__/| \|    \__)| \ \__/|/\| | |  | 
    logstr=sprintf('  Correcting Region Growth (Subj %04d)...\n',isubj);  fprintf(logstr);
    if exist('fid','var'), fprintf(fid,logstr); end
    
    % correct failed slices (0 px in mask) in middle group slices ------------------
    slcs_LVnpx_raw_dia(slcs_LVnpx_raw_dia==0 & midSlcs_msk)=mean(slcs_LVnpx_raw_dia(midSlcs_msk));
    slcs_LVnpx_raw_sys(slcs_LVnpx_raw_sys==0 & midSlcs_msk)=mean(slcs_LVnpx_raw_sys(midSlcs_msk));
    
    % correct basal slices ------------------
    % pick reference slice: (nearest non-basal slice)
    if dirB2A>0,  islice_ref=find(midSlcs_msk,1,'first');
    else          islice_ref=find(midSlcs_msk,1,'last');
    end
    
    % DIA:
    [slcs_LVmsk_dia,slcs_LVnpx_dia]=slcRegGroCorrect(slcs_LVmsk_raw_dia,slcs_LVnpx_raw_dia,...
        islice_ref-dirB2A,-dirB2A,op.areaRatioLim_bas,true,true,slcs_landing_dia);
        % debug view:    [slcs_LVmsk_raw_dia,slcs_LVmsk_dia]; 
        % debug plot:    [slcs_LVnpx_raw_dia,slcs_LVnpx_dia];   
    % SYS:
    [slcs_LVmsk_sys,slcs_LVnpx_sys]=slcRegGroCorrect(slcs_LVmsk_raw_sys,slcs_LVnpx_raw_sys,...
        islice_ref-dirB2A,-dirB2A,op.areaRatioLim_bas,true,true,slcs_landing_sys);
        % debug view:    [slcs_LVmsk_raw_sys,slcs_LVmsk_sys];  
        % debug plot:    [slcs_LVnpx_raw_sys,slcs_LVnpx_sys];   
    
    % correct apical slices ------------------
    % (inputs are basal-corrected slices, NOT raw slices) 
    % pick reference slice: (nearest non-apical slice)
    if dirB2A>0,  islice_ref=find(midSlcs_msk,1,'last');
    else          islice_ref=find(midSlcs_msk,1,'first');
    end
    
    % DIA:
    [slcs_LVmsk_dia,slcs_LVnpx_dia]=slcRegGroCorrect(slcs_LVmsk_dia,slcs_LVnpx_dia,...
        islice_ref+dirB2A, dirB2A, op.areaRatioLim_api, true, false);
        % debug view:    [slcs_LVmsk_raw_dia,slcs_LVmsk_dia];  
        % debug plot:    [slcs_LVnpx_raw_dia,slcs_LVnpx_dia];   
    % SYS:
    [slcs_LVmsk_sys,slcs_LVnpx_sys]=slcRegGroCorrect(slcs_LVmsk_sys,slcs_LVnpx_sys,...
        islice_ref+dirB2A, dirB2A, op.areaRatioLim_api, true, false);
        % debug view:    [slcs_LVmsk_raw_sys,slcs_LVmsk_sys];  
        % debug plot:    [slcs_LVnpx_raw_sys,slcs_LVnpx_sys];   
    
    % dump mp4
    %      __           __          __  
    % |\/||__) |__|    |  \/  \|\/||__) 
    % |  ||       |    |__/\__/|  ||    
    
    % print
    logstr=sprintf('  Saving MP4 (Subj %04d)...\n',isubj);  fprintf(logstr);
    if exist('fid','var'), fprintf(fid,logstr); end
    
    reportMP4RegGro=[reportMP4RegGroRaw,[slcs_LVmsk_dia;slcs_LVmsk_sys]];
    % [dia frame, dia region growth, dia region corrected;
    %  sys frame, sys region growth, sys region corrected]
    dumpMP4_img3D(reportMP4RegGro,sprintf('%s%s%04d regs [dia%02d sys%02d]',saveFolder_top,filesep,isubj,subj_iFr_dia,subj_iFr_sys),10);
    
    %% Clean area data slices w nearly identical locations
    slcs_loc_rnd=round(slcs_loc,4);
    slcs_loc_rnd_u=unique(slcs_loc_rnd);
    nslice_u=length(slcs_loc_rnd_u);
    slcs_LVmm2_dia=slcs_LVnpx_dia.*prod(slcs_pxsize,2);
    slcs_LVmm2_sys=slcs_LVnpx_sys.*prod(slcs_pxsize,2);
    slcs_LVmm2_mono_dia=zeros(nslice_u,1);
    slcs_LVmm2_mono_sys=zeros(nslice_u,1);
    for islice=1:nslice_u
        msk = slcs_loc_rnd==slcs_loc_rnd_u(islice);
        slcs_LVmm2_mono_dia(islice)=mean(slcs_LVmm2_dia(msk));
        slcs_LVmm2_mono_sys(islice)=mean(slcs_LVmm2_sys(msk));
    end
    % debug plot:
    %     figure; plot(slcs_loc,slcs_LVmm2_dia,'s-');
    %     figure; plot(slcs_loc,slcs_LVmm2_sys,'s-');
    %     figure; plot(slcs_loc_rnd_u,slcs_LVmm2_mono_dia,'s-');
    %     figure; plot(slcs_loc_rnd_u,slcs_LVmm2_mono_sys,'s-');
    subj_LVml_dia=trapz(slcs_loc_rnd_u,slcs_LVmm2_mono_dia)/1000;
    subj_LVml_sys=trapz(slcs_loc_rnd_u,slcs_LVmm2_mono_sys)/1000;
    fprintf(fid2,'%d,%.4f,%.4f\n',isubj,subj_LVml_sys,subj_LVml_dia);
    
end % subj loop
fclose(fid); clear fid;
fclose(fid2); clear fid2;



