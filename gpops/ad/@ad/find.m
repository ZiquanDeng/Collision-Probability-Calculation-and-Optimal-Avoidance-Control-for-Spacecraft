function [ix,jx,sx] = find(x, varargin)

% AD implementation of find.m
% Code written by David Benson
% August 2011


if nargout == 1
  if nargin == 1
    ix = find(x.value);
  elseif nargin == 2
    if length(varargin) == 1
      ix = find(x.value, varargin{1});
    else
      ix = find(x.value, varargin{1}, varargin{2});
    end
  end
else
  if nargin == 1
    [ix,jx] = find(x.value);
  elseif nargin == 2
    if length(varargin) == 1
      [ix,jx] = find(x.value, varargin{1});
    else
      [ix,jx] = find(x.value, varargin{1}, varargin{2});
    end
  end
  if nargout == 3
    s.type = '()'; s.subs = {sub2ind(size(x.value),ix,jx)}; sx = subsref(x,s);
  end
end
