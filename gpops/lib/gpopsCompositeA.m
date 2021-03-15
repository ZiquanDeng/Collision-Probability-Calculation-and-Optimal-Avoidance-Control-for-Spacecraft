function A = gpopsCompositeA(Asect)
%------------------------------------------------------------------%
% This function computes a composite Radau pseudospectral          %
% integration matrix.                                              %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%
sections = size(Asect,2);

for i=1:sections
    nodes(i) = size(Asect{i},1);
end

A=sparse(zeros(sum(nodes),sum(nodes)));
rowshift = 0;
colshift = 0;
for i=1:sections
    rowStart  = rowshift+1;
    rowFinish = rowshift+nodes(i);
    columnStart  = colshift+nodes(i);
    columnFinish = colshift+nodes(i);
    rowIndices = rowshift+1:rowshift+nodes(i);
    columnIndices = colshift+1:colshift+nodes(i);
    A(rowIndices,columnIndices) = Asect{i};
    rowshift = rowshift+nodes(i);
    colshift = colshift+nodes(i);
end;
