function D = gpopsCompositeD(Dsect)
%------------------------------------------------------------------%
% This function computes a composite Radau pseudospectral          %
% differentiation matrix.                                          %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%
sections = size(Dsect,2);

for i=1:sections
    nodes(i) = size(Dsect{i},1);
end

D=sparse(sum(nodes),sum(nodes)+1);
rowshift = 0;
colshift = 0;
for i=1:sections
    rowStart  = rowshift+1;
    rowFinish = rowshift+nodes(i);
    columnStart  = colshift+nodes(i);
    columnFinish = colshift+nodes(i)+1;
    rowIndices = rowshift+1:rowshift+nodes(i);
    columnIndices = colshift+1:colshift+nodes(i)+1;
    D(rowIndices,columnIndices) = Dsect{i};
    rowshift = rowshift+nodes(i);
    colshift = colshift+nodes(i);
end;
