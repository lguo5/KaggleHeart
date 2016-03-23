function out=softstep(x,xt,w)
% 
% Generates a soft step function where transition is sigmoidal (sharpness is adjustable).
% 
% INPUT
% 
%     x: vector of x coordinates, the grid on which output is built,
%        e.g. 0:599
%        You can pass in x of multiple different rows, as long as all rows have the same length.
%     
%     xt: Scalar or vector of same length as there are rows in x, in which case one element works on a row in x.
%         In a hard step function, at x=xt the output transitions to 1 from 0.
%         In this soft step output, at x=xt the transition goes to 0.5 (half way).
%         It need not match any value in x, i.e. it can fall between x grid points.
% 
%     w: Transition width: x distance it takes for the output to sigmoidally transition
%        from 0.5 to 0.9933 (to the right) and 0.0067 (to the left). Defaults to 1. 
%        Scalar or vector of same length as there are rows in x, in which case one element works on a row in x.
% 
% OUTPUT
% 
%     out: each row is a soft step function, a vector of same length as x, where at x=xt the output transitions to 0.5.
% 
ncol=size(x,2);
nrow=size(x,1);
if nargin<3 || isempty(w), w=1; end

if nrow>1
    if length(xt)==1, xt=repmat(xt,[nrow 1]); end
    if length(w) ==1,  w=repmat(w, [nrow 1]); end
    if length(xt)==nrow, xt=repmat(xt(:),[1 ncol]); end
    if length(w) ==nrow,  w=repmat(w(:), [1 ncol]); end
end
out=1./(1+exp(-(x-xt)./(w./5))); % /5 so that when w=1 the sigmoid fxn go to 0.9933, 0.0067 in x distance of 1.

% tests: plot one of these:
%     softstep([1:100;1:100],[20;80],3)';
%     softstep([1:100;1:100],[20;80],[3 5])';
%     softstep([linspace(-5,2,100);linspace(-5,2,100)],[-1;-2])';