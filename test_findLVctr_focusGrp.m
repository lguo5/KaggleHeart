
%% Define target list (tg): focus group (hand-picked):
% tg={...
%     '/Users/liheng/Datasets/DSB16/train/1/study/sax_12';... % mid slice, normal
%     '/Users/liheng/Datasets/DSB16/train/1/study/sax_16';... % mid-apical, blurry
%     '/Users/liheng/Datasets/DSB16/train/1/study/sax_6';... % basal, LV open
%     '/Users/liheng/Datasets/DSB16/train/1/study/sax_5';... % heart volves
%     '/Users/liheng/Datasets/DSB16/train/2/study/sax_5';... % heart volves
%     '/Users/liheng/Datasets/DSB16/train/14/study/sax_17';... % hypertrophy, mid-apical
%     '/Users/liheng/Datasets/DSB16/train/14/study/sax_21';... % hypertrophy, most apical, 
%     '/Users/liheng/Datasets/DSB16/train/3/study/sax_50';... % apical, very small
%     };
% saveFolder='focusGroup';

%% Define target list (tg): all mid-slice cines

% path_group='/Users/liheng/Datasets/DSB16/train';
% subjnums=1:500;
% saveFolder='midSlcTest_train';

% path_group='/Users/liheng/Datasets/DSB16/validate';
% subjnums=501:700;
% saveFolder='midSlcTest_validate';

path_group='/Users/liheng/Datasets/DSB16/test';
subjnums=701:1140;
saveFolder='midSlcTest_test';

nsubj=length(subjnums);
tg=cell(nsubj,1);
for isubj=1:nsubj
    path_study=[path_group,filesep,num2str(subjnums(isubj)),filesep,'study',filesep];
    list_cines=dir([path_study,'sax*']);
    ncine=length(list_cines);
    cineSliceLoc=nan(1,ncine);
    for icine=1:ncine
        path_cine=[path_study,list_cines(icine).name,filesep];
        list_dcm=dir([path_cine,'*.dcm']);
        ndcm=length(list_dcm);
        if ndcm<40
%             warn('Subj %d -> "%s" contains %d DICOM frames. Subject will be skipped...',subjnums(isubj),pathLastPart(path_cine),ndcm);
%         else
            [~, dcmInfo]=loadCine(path_cine,'frames',1,'print',false,'dcmInfoOnly',true);
            cineSliceLoc(icine)=dcmInfo{1}.SliceLocation;
        end
    end
    % skips subj who has any cine with too many frames
    if any(isnan(cineSliceLoc))
        indMidSlc=nan;
        tg{isubj}='';
        fprintf('Subj %4d:  SKIPPED. slcloc=%s\n',subjnums(isubj),mat2str(cineSliceLoc,4));
    else
        [~,indMidSlc]=min(abs(cineSliceLoc-median(cineSliceLoc)));
        tg{isubj}=[path_study,list_cines(indMidSlc).name];
        fprintf('Subj %4d:  MidSlc=%2d="%s"\tslcloc=%s\n',subjnums(isubj),indMidSlc,list_cines(indMidSlc).name,mat2str(cineSliceLoc,4));
    end
end

%% load
% ncine=length(tg);
% im3d_raw_all=cell(ncine,1);
% for icine=1:ncine
%     % load cine
%     fprintf('Loading %s...\n',tg{icine});
%     [im3d_raw_all{icine}, dcmInfo]=loadCine(tg{icine},'print',false);    
% end

%% start
ncine=length(tg);
if ~isempty(saveFolder),mkdir(saveFolder); end

for icine=1:ncine
    
    if ~isempty(tg{icine});
        % find subj and cine names:
        [ia,ib]=regexp(tg{icine},'/\d*/');
        name_subj=tg{icine}(ia+1:ib-1);
        [ia,ib]=regexp(tg{icine},'sax_\d*');
        name_cine=tg{icine}(ia+4:ib);
        fprintf('Finding LV center: Subj %s, sax_%s...\n',name_subj,name_cine);

        [im3d_raw, dcmInfo]=loadCine(tg{icine},'print',false);    
        % im3d_raw=im3d_raw_all{icine};

        % crop -------
        n1=size(im3d_raw,1);
        n2=size(im3d_raw,2);
        im3d=im3d_raw(round(n1/4)+1:round(n1/4*3),round(n2/4)+1:round(n2/4*3),:);
        % im3d=im3d_raw;

        im3d=LG_normalize(im3d,[],[0 99]);
        n1=size(im3d,1);
        n2=size(im3d,2);
        n3=size(im3d,3);

        %% separate myocardium, blood:
        [msk3d_myo,msk3d_bld, intens_myo,intens_bld, msk3d_hrt,im3d_mopx]=cineClusterBldMyo(im3d);
        % clean up a bit:
        msk3d_myo=runFxnOnPage(msk3d_myo,@medfilt2, [3 3]); %%%%
        msk3d_bld=runFxnOnPage(msk3d_bld,@medfilt2, [3 3]); %%%%
        msk3d_myo=runFxnOnPage(msk3d_myo,@maskTrimRegions, 0.9, 4); %%%%
        msk3d_bld=runFxnOnPage(msk3d_bld,@maskTrimRegions, 0.9, 4); %%%%
    %     msk3d_bld_ero=runFxnOnPage(msk3d_bld,@imerode,strel('disk',2));
    %     msk3d_bld_ero=runFxnOnPage(msk3d_bld_ero,@imdilate,strel('disk',2));
    %     runFxnOnPage(msk3d_myo,@centroidOfMask)
    %     runFxnOnPage(msk3d_myo,@maskTrimRegions, 0.9, 4);
    %     runFxnOnPage(msk3d_myo,@medfilt2, [3 3]);

        %% blood region roundness 
        rMap3d_bld=runFxnOnPage(msk3d_bld,@maskRegRndmap,'convex',true);  %   sum(rMap3d_bld,3); 
        rMap3d_bld2=rMap3d_bld; % sum(rMap3d_bld,3); 
        rMap3d_bld2(rMap3d_bld2<0.9*max(rMap3d_bld(:)))=0; %%%%   %   sum(rMap3d_bld2,3); 

        %% bwdist on inverse blood mask
        bwdist3d=runFxnOnPage(~msk3d_bld,@bwdist);
        bwdist3d=LG_scale(bwdist3d);
        
        %% Mark roundness map max
        rMap2d=LG_scale(mean(rMap3d_bld2,3));
        c1c2=matrixFindMaxRegCtr(rMap2d,0.9); %%%%
        rMap2d(c1c2(1),:)=1;
        rMap2d(:,c1c2(2))=1;
        
        %% mixed BWdist and roundness
        RnBW2d=LG_scale(mean(bwdist3d.*rMap3d_bld2,3));
        c1c2=matrixFindMaxRegCtr(RnBW2d,0.9); %%%%
        RnBW2d(c1c2(1),:)=1;
        RnBW2d(:,c1c2(2))=1;
        
        %% watershed
    %     [D,IDX] = bwdist(~msk3d_bld(:,:,1));  % -D  
    %     D=-D;
    %     D(~msk3d_bld(:,:,1))=-Inf;
    %     L = watershed(D);


    %     %% run rays 
    %     acum=zeros(size(im3d));
    %     gmag=zeros(size(im3d));
    %     gdir=zeros(size(im3d));
    %     [a1, a2]=ind2sub([n1,n2],1:(n1*n2));
    % 
    %     im4grad=convn(msk3d_bld,ones(3,3,3),'same');
    %     [gx,gy,gz]=gradient(im4grad); %    gx.^2+gy.^2+gz.^2;   abs(gz);   
    %     % (gx.^2+gy.^2).*gz.^2;    %   gmag.*gz.^2;
    %     
    %     for ifr=1:n3
    %         [gmag(:,:,ifr),gdir(:,:,ifr)]=imgradient(im4grad(:,:,ifr));  %  
    %     end
    %     for ifr=1:n3
    %         gmag_curr=gmag(:,:,ifr);
    %         gdir_curr=gdir(:,:,ifr);
    %         msk=gmag_curr>prctile(gmag_curr(:),80); %     msk=msk_gmag(:,:,ifr);
    %         acum(:,:,ifr)=drawRay(zeros(n1,n2),a1(msk(:)),a2(msk(:)),gdir_curr(msk(:)),gmag_curr(msk(:)),30,true);
    %     end


        % voting pool: rays, bwdist, region roundness, polar image 

%         %% save mp4
%         im_cat=[im3d,im3d_mopx,msk3d_myo;msk3d_bld,repmat(rMap2d,[1 1 n3]),repmat(RnBW2d,[1 1 n3])];
%         dumpMP4_img3D(im_cat, sprintf('%s%s%04d sax_%s',saveFolder,filesep,str2double(name_subj),name_cine), 30, 80)
    
        %% save image
        im_cat=[im3d(:,:,1),im3d_mopx(:,:,1),msk3d_myo(:,:,1);msk3d_bld(:,:,1),rMap2d,RnBW2d];
        imwrite(im_cat,sprintf('%s%s%04d sax_%s.jpg',saveFolder,filesep,str2double(name_subj),name_cine),'jpg','quality',80);
        
    end % if current tg is not empty
    
end % looping thru cines





%% junk

    % msk3d_bld_wk=convn(msk3d_bld,ones(3,3,1),'same'); %    abs(diff(msk3d_bld_wk,1,3));   

%     %% vector containment
%     gConvKer=ones(24,24,1)./24^2;
%     gConvBorder='valid';
%     gsum=convn(gx,gConvKer,gConvBorder).^2+convn(gy,gConvKer,gConvBorder).^2;
%     gmag2=convn(gmag,gConvKer,gConvBorder);
%     numer=LG_scale(gmag2);
%     denom=LG_scale(gsum);
%     numer(numer<=0)=eps;
%     denom(denom<=0)=eps;
%     numer-denom



