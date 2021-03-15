function y = compositeDerivative(x,y,outerDerivative)

nDerivatives = x.nderivs;
leadingDim = numel(x.value);
[rows,cols,innerDerivative] = find(x.derivative);
derivative = outerDerivative(rows).*innerDerivative;
y.derivative = sparse(rows,cols,derivative,leadingDim,nDerivatives);
y.nderivs = nDerivatives;
