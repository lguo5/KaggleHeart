function out=runFxnOnPage(img3d,fhandle,varargin)
% Runs a function (fhandle) on each page of a 3D array, passing varargin to the function fhandle.
% e.g. im3d_out=runFxnOnPage(im3d_in, @medfilt2, [3 3]);

% out=img3d; % remain type agnostic
n3=size(img3d,3);

% test apply:
out1=fhandle(img3d(:,:,1),varargin{:});
out=repmat(out1,[1 1 n3]);

for ifr=2:n3
    out(:,:,ifr)=fhandle(img3d(:,:,ifr),varargin{:});
end

% All outputs of fh will be concat'd together along the next higher
% dim of the single output. e.g. if a single output is [a b] or [a;b],
% out1=fh(img3d(:,:,1),varargin{:});
% nd=ndims(out1);
% out=repmat(out1,[ones(1,nd),n3]);
