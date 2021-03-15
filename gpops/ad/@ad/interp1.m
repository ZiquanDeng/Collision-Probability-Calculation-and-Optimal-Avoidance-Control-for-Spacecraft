function yi = interp1(x,y,xi_ad,varargin)

% AD implementation of interp1.m
% Code written by David Benson
% Aug 2011
%
% Note:  This function currently only allows 'spline', 'cubic', 
%        'phchip', or 'linear' interpolation.  The
%        other methods will produce an error
%       (see documentation for MATLAB function INTERP1)
% Inputs and Outputs
%      x: known vector of x's for interpolation
%      y: known vector of y's for interpolation
%      xi: object of AD class 
%      method: 'spline', 'cubic', 'phchip', or 'linear'
%      yi: object of AD class 

if (nargin<3) || ...
   (nargin==3 && ischar(xi_ad))  || ...
   (nargin==4 && ~(ischar(varargin{1}) || isempty(varargin{1})) || ...
   (nargin==4 && strcmp(varargin{1}, 'extrap'))) || ...
   (nargin>4);
    error('GPOPS:interp1:InvalidCall','interp1 must be called YI = INTERP1(X,Y,XI,METHOD) for automatic derivative')
end

%% process Y
siz_y = size(y);
% y may be an ND array, but collapse it down to a 2D yMat. If yMat is
% a vector, it is a column vector.
if isvector(y)
    if isrow(y)
        % Prefer column vectors for y
        yMat = y.';
        n = siz_y(2);
    else
        yMat = y;
        n = siz_y(1);
    end
    ds = 1;
    prodDs = 1;
else
    error('GPOPS:interp1:Yvector','Y must be a vector for automatic derivative.');
%     n = siz_y(1);
%     ds = siz_y(2:end);
%     prodDs = prod(ds);
%     yMat = reshape(y,[n prodDs]);
end

%% process X
if ~isvector(x)
    error('MATLAB:interp1:Xvector','X must be a vector.');
end
if length(x) ~= n
    if isvector(y)
        error('MATLAB:interp1:YInvalidNumRows', ...
            'X and Y must be of the same length.')
    else
        error('MATLAB:interp1:YInvalidNumRows', ...
            'LENGTH(X) and SIZE(Y,1) must be the same.');
    end
end
% Prefer column vectors for x
xCol = x(:);

%% Process XI
xi = xi_ad.value;
siz_xi = size(xi);
% xi may be an ND array, but flatten it to a column vector xiCol
xiCol = xi(:);
% The size of the output YI
if isvector(y)
    % Y is a vector so size(YI) == size(XI)
    siz_yi = siz_xi;
else
    if isvector(xi)
        % Y is not a vector but XI is
        siz_yi = [length(xi) ds];
    else
        % Both Y and XI are non-vectors
        siz_yi = [siz_xi ds];
    end
end

%% check inputs
if ~isreal(x)
    error('MATLAB:interp1:ComplexX','X should be a real vector.')
end

if ~isreal(xi)
    error('MATLAB:interp1:ComplexInterpPts', ...
        'The interpolation points XI should be real.')
end

% Error check for NaN values in X and Y
% check for NaN's
if (any(isnan(xCol)))
    error('MATLAB:interp1:NaNinX','NaN is not an appropriate value for X.');
end

% NANS are allowed as a value for F(X), since a function may be undefined
% for a given value.
if any(isnan(yMat(:)))
    warning('MATLAB:interp1:NaNinY', ...
        ['NaN found in Y, interpolation at undefined values \n\t',...
        ' will result in undefined values.']);
end

if (n < 2)
    error('MATLAB:interp1:NotEnoughPts', ...
        'There should be at least two data points.')
end

%% Process METHOD in
if nargin >= 4 && ~isempty(varargin{1})
    method = varargin{1};
else
    method = 'linear';
end

%% start algorithm

h = diff(xCol);
if any(h < 0)
    [xCol,p] = sort(xCol);
    yMat = yMat(p,:);
    h = diff(xCol);
end
if any(h == 0)
    error('MATLAB:interp1:RepeatedValuesX', ...
        'The values of X should be distinct.');
end

% Interpolate
numelXi = length(xiCol);
switch method(1)
    case 's'  % 'spline'
        % get spline poly values
        pp = spline(xCol,yMat);

        %eval values
        yiMat = ppval(pp,xiCol);

        % compute derivative poly values
        ppDer = pp;
        numCoef = size(pp.coefs,2);
        ppDer.coefs(:,1) = 0;
        for i = 1:numCoef-1
           ppDer.coefs(:,i+1) = (numCoef - i) * pp.coefs(:,i);
        end

        % eval derivatives
        outerDerivative = ppval(ppDer,xiCol);

    case {'c','p'}  % 'cubic' or 'pchip'
        % get poly values
         pp = pchip(xCol,yMat);

        %eval values
        yiMat = ppval(pp,xiCol);

        % compute derivative poly values
        ppDer = pp;
        numCoef = size(pp.coefs,2);
        ppDer.coefs(:,1) = 0;
        for i = 1:numCoef-1
           ppDer.coefs(:,i+1) = (numCoef - i) * pp.coefs(:,i);
        end

        % eval derivatives
        outerDerivative = ppval(ppDer,xiCol);
        
  case{'l'} % linear
        yiMat = zeros(numelXi,prodDs,superiorfloat(xCol,yMat,xiCol));
        outerDerivative = zeros(numelXi,prodDs,superiorfloat(xCol,yMat,xiCol));
        
        p = 1:numelXi;
        
        [~,k] = histc(xiCol,xCol);
        k(xiCol<xCol(1) | ~isfinite(xiCol)) = 1;
        k(xiCol>=xCol(n)) = n-1;
        
        s = (xiCol - xCol(k))./h(k);
        
        for j = 1:prodDs
            yiMat(p,j) = yMat(k,j) + s.*(yMat(k+1,j)-yMat(k,j));
            outerDerivative(p,j) = (yMat(k+1,j)-yMat(k,j))./h(k);
        end
               
  otherwise % 'nearest','v5cubic'
     error('GPOPS:interp1:InvalidMethod',...
       ['Invalid method. Only ''spline'', ''cubic'', ''phchip'', ',...
       'or ''linear'' interpolation supported for automatic derivatives']) 
end

yi.value = reshape(yiMat,siz_yi);
yi = compositeDerivative(xi_ad,yi,outerDerivative(:));
yi = class(yi,'ad');

