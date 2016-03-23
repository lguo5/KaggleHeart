function [fileName, filePath] = pathLastPart(fullPath)
% Returns the trailing folder name or the trailing file name (Returns
% everything after the last filesep. If incoming path ends with a filesep,
% that filesep is stripped first.)
% 
% '/some/folders/foo'  --> 'foo', '/some/folders'
% '/some/folders/foo/' --> 'foo', '/some/folders'
% '/some/folders/foo/001.txt' --> '001.txt', '/some/folders/foo'
% 

if strcmpi(fullPath(end),filesep), fullPath=fullPath(end-1); end % if the last char is a filesep, strip it.

lastFileSepInd = strfind( fullPath, filesep);
lastFileSepInd = lastFileSepInd(end);

fileName = fullPath( lastFileSepInd+1 : end );
filePath = fullPath( 1:lastFileSepInd-1 ); 