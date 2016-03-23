function [im3d, dcmInfo]=loadCine(targetpath,varargin)
% Loads one cine, into a 3D image.
% 
%     [im3d, dcmInfo]=loadCine(targetpath,'frames',3:5,'print',true,'dcmInfoOnly',false)
% 
% INPUT
%     
%     targetpath: e.g. '~/Dropbox/KaggleBowl/1/study/sax_5'
%     
% INPUT: PARAMETER-VALUE PAIRS:
%     
%     frames: list of dicom frames to load, e.g. 1, [1,30], 1:30... Use [] to load all.
% 
%     dcmInfoOnly: (Defaults to False) If TRUE only loads dicom tag data, returns all-0 for actual image. Default=False.
% 
%     dcmInfo1stFrOnly: (Defaults to F) If T will load dicom tags for 1st image.
% 
%     print: if TRUE (default), prints which frame is being loaded.
%         
% OUTPUT
% 
%     dcmInfo:
%         Cell array of DICOM tags, each element (a struct) is for a DICOM image.
%         If dcmInfo1stFrOnly is True, dcmInfo will be a one-element cell.
%     
%     

[frames,varargin]=findStripArg(varargin,'frames',[]);
[dcmInfoOnly,varargin]=findStripArg(varargin,'dcmInfoOnly',false);
[dcmInfo1stFrOnly,varargin]=findStripArg(varargin,'dcmInfo1stFrOnly',false);
[toprint,varargin]=findStripArg(varargin,'print',true);

dcmList = dir([targetpath filesep '*.dcm']);
dcmInfo1 = dicominfo([targetpath filesep dcmList(1).name]);

% sizes for dim1~3
% n1=dcmInfo.Height;
% n2=dcmInfo.Width;
n1=dcmInfo1.Rows;
n2=dcmInfo1.Columns;
n3=length(dcmList);
if isempty(frames), 
    frames=1:n3;
else
    n3=length(frames);
end

im3d=zeros(n1,n2,n3);

% dcmInfo=repmat(dcmInfo1,[n3,1]);
if dcmInfo1stFrOnly,
    dcmInfo=cell(1,1);
else
    dcmInfo=cell(n3,1);
end

if toprint, fprintf('Loading %d DICOMs from ...%s, frame: ',n3,targetpath(end-9:end)); end
for ind=1:n3
    if toprint, fprintf('%d ',frames(ind)); end
    if ~dcmInfoOnly
        im3d(:,:,ind)=dicomread([targetpath filesep dcmList(frames(ind)).name]);
    end
    if ~dcmInfo1stFrOnly || (dcmInfo1stFrOnly && ind==1)
        dcmInfo{ind} =dicominfo([targetpath filesep dcmList(frames(ind)).name]);
    end
end
if toprint, fprintf('Done.\n'); end


% Test: 
%      [im3d, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_12');    
%      [im3d, dcmInfo]=loadCine('/Users/liheng/Dropbox/KaggleBowl/1/study/sax_12','frames',1,'print',false);    
% 
