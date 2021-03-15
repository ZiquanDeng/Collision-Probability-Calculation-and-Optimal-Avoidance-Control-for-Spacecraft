function y = prod(x,dim)

% AD implementation of prod.m
% Code written by David Benson
% August 2011

n = size(x.value);
if nargin == 1  
  dim = find(n > 1, 1, 'first');
  if isempty(dim)
    y = x;
    return
  end
end

if dim > length(n)
  error('GPOPS:prod:InvalidDim','invalid dimension input')
end

nnew = n;
nnew(dim) = 1;
y = ones(nnew);

s.type = '()'; s.subs = cell(1,ndims(x.value));
for i = 1:ndims(x.value)
  s.subs{i} = ':';
end
for i = 1:n(dim)
  s.subs{dim} = i; 
  xnew = subsref(x,s);
  
  y = y .* xnew; 

end

