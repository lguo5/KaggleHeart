function out=getAgeYr(dcmInfo)
% Takes in a structure containing dicom tag info as returned by Matlab
% dicominfo(), finds the patient age tag ('000Y','000M',etc), returns in num of years.

ageStr=dcmInfo.PatientAge;
if isempty(regexp(ageStr,'\d\d\d\w', 'once')),
    error('%s: Unrecognized age format.',mfilename);
end
out=str2double(ageStr(1:3));
if strcmpi(ageStr(end),'M')
    out=out/12;
elseif strcmpi(ageStr(end),'W')
    out=out*7/365;
elseif strcmpi(ageStr(end),'D')
    out=out/365;
end
