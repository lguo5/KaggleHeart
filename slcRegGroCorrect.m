function [mskIn,areaIn]=slcRegGroCorrect(mskIn,areaIn,iCurr,indIncr,areaRatioLim,fillZero,doWatershed,ctrs)
%     mskIn: stacks of region masks.
%     ctrs: [nfr 2], region grow seed point for each slice, each row is [d1 d2].
%     areaRatioLim: e.g. 1.1, ratio of CURRENT slice region area to PREV slice area.
%     indIncr: +1 or -1, meaning to go forward or backward along dim3 of input 3D mask.
%     iCurr: starting TARGET slice ind, not reference slice ind
n3=size(mskIn,3);

% mskOut=mskIn;
% areaOut=areaIn;

while iCurr>=1 && iCurr<=n3
    % disp(iCurr); % debug print
    iPrev=iCurr-indIncr;
    
    % correct 0 if enabled:
    if fillZero && areaIn(iCurr)==0
        mskIn(:,:,iCurr)=mskIn(:,:,iPrev);
        areaIn(iCurr)=areaIn(iPrev); 
    
    % do watershed if enabled:
    elseif areaIn(iCurr)/areaIn(iPrev)>areaRatioLim
        if doWatershed
            D=-bwdist(~mskIn(:,:,iCurr));
            D(~mskIn(:,:,iCurr))=-inf;
            L=watershed(D);
            stats = regionprops(L,{'Area','Centroid','Perimeter'});
            
            % only keep those area with more than 10 px
            stats=stats(cell2mat({stats(:).Area})>10);
            nreg=length(stats);
            
            % Choose region with min distance to original seed point:
%             reg_dist2ctr=sum(( cell2mat({stats(:).Centroid}') - repmat(ctrs(iCurr,:),[nreg,1]) ).^2,2);
%             [~, isort1]=sort(reg_dist2ctr); % isort: small=better
            
            % Choose region w largest roundness
            reg_area=cell2mat({stats(:).Area});
            reg_peri=cell2mat({stats(:).Perimeter});
            reg_roundness=4.*pi.*reg_area/reg_peri.^2;
            [~, isort2]=sort(reg_roundness,'descend'); % isort: small=better
            % [~, imin]=max(reg_roundness);
            
            % combine
            isort=isort2;
            [~, imin]=min(isort);
            
            mskIn(:,:,iCurr)=L==imin;
            areaIn(iCurr)=stats(imin).Area;
        end
        if areaIn(iCurr)/areaIn(iPrev)>areaRatioLim % if condition persists, just wrap:
            mskIn(:,:,iCurr)=mskIn(:,:,iPrev);
            areaIn(iCurr)=areaIn(iPrev); 
        end
    end
    iCurr=iCurr+indIncr;
end