%% Intro
% Loops through all SAX folders of all cases in a group to check/print some
% property of each cine.
% 
% HOW TO USE:
% >> To read from DICOMs, run tests and save .mat of all DICOM info, set
% fromDicom=TRUE.
% >> To skip the DICOMs and run tests off of the saved .mat (much faster), set
% fromDicom=False. Expects user having dragged an existing .mat to
% workspace, containing: all_dcmInfo, all_cinedirs, all_subjdirs

%% select target folder (of all subjects)
fromDicom=false;

frames=1; % when reading from dicoms, num of frames to load for each cine. Use [] to load all. Sometimes you just need tags from 1st frame.

% path_group='/Users/liheng/Datasets/DSB16/train';
% path_group='/Users/liheng/Datasets/DSB16/validate';
path_group='/Users/liheng/Datasets/DSB16/test';

% saveName='Image Size Greater than 256.txt'; % use '' to skip text file saving
% saveName='All Image Size Width x Height.txt'; % use '' to skip text file saving
% saveName='RowNotEqualHei or ColNotEqualWid.txt'; % use '' to skip text file saving
% saveName='Image Size Unequal.txt';
saveName='Pixel Spacing Unequal.txt';
% saveName='Image Thickness Unequal.txt';
% saveName='All SAX folders.txt';
% saveName='nPhase not 30.txt';
% saveName='All Patient Sex Strings.txt';
% saveName='LUT Folder Name vs Nth Loaded.txt';


%% start
if fromDicom
    all_subjdirs=dir(path_group);
    % remove anything whose name starts with . or isn't a directory
    maskExclude=cellfun(@(x)strcmp(x(1),'.'),{all_subjdirs.name}) | cellfun(@(x)~x,{all_subjdirs.isdir}); % first define an anonymous fxn that compares 1st char of input to '.', then apply it to all dir names forced into a cell.
    all_subjdirs=all_subjdirs(~maskExclude);
end
nsubj=length(all_subjdirs);


if fromDicom
    % {subj1;subj2;...}
    % where subjX={DicomInfoStructArrayCine1;DicomInfoStructArrayCine2;...},
    % where Nth subj corresponds to Nth in all_subjdirs
    all_dcmInfo=cell(nsubj,1); 

    % {dirStructArraySubj1;dirStructArraySubj2;...}, where
    % subj are in same order as in all_subjdirs.
    all_cinedirs=cell(nsubj,1); 
end

writeText=~isempty(saveName);
if writeText, fid=fopen(saveName,'w');  end

% % For slice location distribution plot ------------
% figure_LG('Slice Locations',[20 200 1200 480]);

for isubj=1:nsubj
    fprintf('Subject folder: %s\n',all_subjdirs(isubj).name);
    
    if fromDicom
        path_subj=[path_group,filesep,all_subjdirs(isubj).name,filesep,'study',filesep];
        list_cine=dir([path_subj,'sax*']); % <<<< SAX filter
    else
        list_cine=all_cinedirs{isubj};        
    end    
    
    ncine=length(list_cine);
    if ncine==0
        warning('!!!! %s does not have STUDY folder\n', all_subjdirs(isubj).name);
        % all_dcmInfo{isubj} will be left empty
    else
        if writeText
            fprintf(fid,'%s: ',all_subjdirs(isubj).name); % prints beginning of a line: 1: 2: ...
        end
        if fromDicom
            all_cinedirs{isubj}=list_cine; % record sax folder list
            dcmInfo_1subj=cell(ncine,1);
        end
        
%         % For 'Image Size Unequal.txt'; -----------------
%         subj_w=zeros(ncine,1);
%         subj_h=zeros(ncine,1);

        % For 'Pixel Spacing Unequal.txt'; -----------------
        subj_pxspa=zeros(ncine,2);

%         subj_t=zeros(ncine,1);
%         subj_z=zeros(ncine,1); % For slice location distribution plot ------------
        
%         % For 'All Patient Sex Strings.txt'; -----------------
%         if writeText,fprintf(fid,'%s',all_dcmInfo{isubj}{1}.PatientSex); end

%         % For 'LUT Folder Name vs Nth Loaded.txt'; -----------------
%         if writeText,fprintf(fid,'%d',isubj); end

        %% LOOPING THROUGH CINES 
        for icine=1:ncine % sax_1, sax_2, ...
            if fromDicom
                path_cine=[path_subj,list_cine(icine).name,filesep];
                
                [~,dcmInfo]=loadCine(path_cine,'frames',frames,'print',false,'dcmInfoOnly',true);
                dcmInfo_1subj{icine}=dcmInfo; % dcmInfo is struct array when nframe loaded >1
            else
                dcmInfo=all_dcmInfo{isubj}{icine};
            end
            
             
%             % For 'Image Size Greater than 256.txt' --------------
%             maxsize=max([dcmInfo.Rows, dcmInfo.Columns, dcmInfo.Width, dcmInfo.Height]); % maxsize=max(size(im3d));
%             if writeText && maxsize>256, fprintf(fid,'%s=%d ',list_cine(icine).name,maxsize); end
            
%             % Test: if any Rows!=Height, Cols!=Height ---------------
%             if writeText && (dcmInfo.Rows~=dcmInfo.Height || dcmInfo.Columns~=dcmInfo.Width)
%                 fprintf(fid,'%s:Ro/Hi/Co/Wi=%d/%d/%d/%d  ',list_cine(icine).name, dcmInfo.Rows, dcmInfo.Height, dcmInfo.Columns, dcmInfo.Width); 
%             end

%             % For 'All Image Size Width x Height.txt' --------------
%             if writeText,fprintf(fid,'%s:%dx%d  ',list_cine(icine).name, dcmInfo.Width, dcmInfo.Height); end 

%             % For 'Image Size Unequal.txt'; -----------------
%             subj_w(icine)=dcmInfo.Width;
%             subj_h(icine)=dcmInfo.Height;
            
            % For 'Pixel Spacing Unequal.txt'; -----------------
            subj_pxspa(icine,:)=dcmInfo.PixelSpacing(:)';

%             % For 'Image Thickness Unequal.txt'; -----------------
%             subj_t(icine)=dcmInfo.SliceThickness;
            
%             % For 'All SAX folders.txt'; -----------------
%             if writeText,fprintf(fid,'%s  ',list_cine(icine).name); end
            
%             % For 'nPhase not 30.txt'; -----------------
%             if writeText && dcmInfo.CardiacNumberOfImages~=30, fprintf(fid,'%s=%d  ', list_cine(icine).name, dcmInfo.CardiacNumberOfImages); end

%             % For slice location distribution plot -----------------
%             subj_z(icine)=dcmInfo.SliceLocation;

            
        end % cines loop
        
%         % For 'Image Size Unequal.txt'; -----------------
%         if sum(diff(subj_w))>0 || sum(diff(subj_w))>0
%             if writeText, fprintf(fid,'Width=[%s]  Height=[%s]',num2str(subj_w'),num2str(subj_h')); end 
%         else
%             if writeText, fprintf(fid,'All Width x Height = %dx%d  ',dcmInfo.Width, dcmInfo.Height); end 
%         end

        % For 'Pixel Spacing Unequal.txt'; -----------------
        if sum(sum(diff(subj_pxspa,1,1),1),2)>0
            if writeText, fprintf(fid,'PxSpacing(1,2) = %s, %s',mat2str(round(subj_pxspa(:,1)',3)),mat2str(round(subj_pxspa(:,2)',3))); end 
        else
            if writeText, fprintf(fid,'All Pixel Spacing = %s  ',mat2str(round(dcmInfo.PixelSpacing,3))); end 
        end

%         % For 'Image Thickness Unequal.txt'; -----------------
%         if sum(diff(subj_t))>0
%             if writeText, fprintf(fid,'Thickness=[%s]',num2str(subj_t')); end 
%         else
%             if writeText, fprintf(fid,'All Thickness = %d',dcmInfo.SliceThickness); end 
%         end
        
%         % For slice location distribution -----------------
%         plot(ones(ncine,1).*str2double(all_subjdirs(isubj).name),subj_z,'s-'); hold on;
        
        if fromDicom, all_dcmInfo{isubj}=dcmInfo_1subj; end
        if writeText, fprintf(fid,'\n'); end

    end
end
if writeText, fclose(fid); end

% % For slice location distribution -----------------
% xlabel(sprintf('Subjects (%s)',path_group));
% ylabel('Slice Location (mm)');
% zoom xon;

%% save .mat
if fromDicom
    save(sprintf('set=XXXX N=%d frame=XX',nsubj),'all_dcmInfo','all_cinedirs','all_subjdirs');
end
