function mskout=maskTrimRegions(mskin,frxn,nConnect)
% Takes in a mask of various regions, delete smallest regions to keep at least frxn
% fraction of total 1's. useful for deleting small isolated specs in a binary image.
% 
%     nconnect: passed to bwconncomp
%     
regions=bwconncomp(mskin,nConnect);
regions_npx=cellfun(@(x)length(x),regions.PixelIdxList);
[regions_npx_sorted,sortInd]=sort(regions_npx,'descend');
npx_cumsum=cumsum(regions_npx_sorted)./sum(regions_npx_sorted);
iEnough=find(npx_cumsum>=frxn,1,'first');
mskout=false(size(mskin));
if ~isempty(iEnough)
    mskout(cell2mat(transpose(regions.PixelIdxList(sortInd(1:iEnough)))))=true;
end


