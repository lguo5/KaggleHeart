function writeCSV(CDFs,IDs,saveName,formatStr)
% INPUT
% 
%     CDFs: each row (e.g. 600-element) is a CDF.
%         Every 2 rows are for a case: diast and syst rows interleave, e.g.:
%         001_Diastole,...
%         001_Systole,...
%         002_Diastole,...
%         002_Systole,...
% 
%     IDs: vec of case numbers. IDs must have length half as there are rows in CDFs.
%         e.g. [1 2 3...] will be written as:
%         001_Diastole
%         001_Systole
%         002_Diastole
%         002_Systole
%         003_Diastole
%         003_Systole
%         ...
% 
%     saveName: string for the output filename. Include .csv .txt etc. Defaults to SubmitMe.csv
% 
%     formatStr: how each value of CDF is written in: defaults to %.4f.
% 
if nargin<3 || isempty(saveName), saveName='SubmitMe.csv'; end
if nargin<4 || isempty(formatStr), formatStr='%.4f'; end

% nrow=size(CDFs,1);
ncase=length(IDs);
ncol=size(CDFs,2);

fid=fopen(saveName,'w');

% write 1st line: "Id,P0,P1,P2,P3,...,P599
fprintf(fid,'Id');
fprintf(fid,',P%d',0:ncol-1);
fprintf(fid,'\n');

% write 2 rows for each case:
for icase=1:ncase
    % "000_Diastole,0.001,0.002,..."
    fprintf(fid,'%03d_Diastole',IDs(icase));
    fprintf(fid,[',',formatStr],CDFs(icase*2-1,:));
    fprintf(fid,'\n');
    % "000_Systole,0.001,0.002,..."
    fprintf(fid,'%03d_Systole',IDs(icase));
    fprintf(fid,[',',formatStr],CDFs(icase*2,:));
    fprintf(fid,'\n');
    
end

fclose(fid);

