function str=stripNonEnding0(str)
% Delete 0s in a string except 0s that are at the end of the string.
% e.g. 'SC-HF-000I-01000' --> 'SC-HF-I-1000'
inds0=strfind(str,'0');
n=length(str);
while ~isempty(inds0) && inds0(end)==n
    inds0(end)=[];
    n=n-1;
end
str(inds0)=[];

