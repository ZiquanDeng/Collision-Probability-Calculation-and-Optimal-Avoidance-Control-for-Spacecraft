function [xinterp,lagrange] = gpopsLagrange(tinterp,t,x)
%------------------------------------------------------------------%
% This function is for use with the mesh refinement algorithm.     %
% This function calculates arbitrary points on the Lagrange        %
% polynomial defined by the support points at the initial time,    %
% the Legendre-Gauss points and the final time of a segment.       %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%
lagrange = ones(length(tinterp),length(t));
for i = 1:length(t)
    for j = 1:length(t)
        if i ~= j
            lagrange(:,i)=lagrange(:,i).*(tinterp-t(j))./(t(i)-t(j));
        end
    end
end
xinterp = lagrange*x;

% 
% lagrange2 = ones(length(tinterp,length(t));
% 
% for i = 1:length(t)
%     index = find( t ~= t(i) );
%     lagrange2(:,index)=lagrange2(:,index).*(tinterp-i)./(t(index)-i);
% end
