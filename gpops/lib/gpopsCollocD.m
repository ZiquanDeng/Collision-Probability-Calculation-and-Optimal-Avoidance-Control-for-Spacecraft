function [D,Dd,Do]=collocD(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% collocD.m
%
% Computes the pseudospectral/collocation differentiation matrix for the 
% arbitrary nodes stored in the vector x. Uses the lagrange polynomial 
% formulation.
%
% Reference:
%
% Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange Interpolation" 
% http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/berrut.ps.gz
%
% Written by: Greg von Winckel       07/18/04
% Contact:    gregvw@chtm.unm.edu 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make x a column vector if it isn't already and order it
% and get the number of nodes
%x=sort(x(:));                       N=length(x); N1=N+1; N2=N*N;

% Compute the barycentric weights
%X=repmat(x,1,N);                    Xdiff=X-X'+eye(N);
%W=repmat(1./prod(Xdiff,2),1,N);     D=W./(W'.*Xdiff); 
%D(1:N1:N2)=1-sum(D);                D=-D';

%---- Compute Diff matrix

N1 = length(x)-1;
%  xxPlusEnd = [x; 1];
xxPlusEnd = x;
M = length(xxPlusEnd);
M1 = M+1;
M2 = M*M;

% Compute the barycentric weights
Y = repmat(xxPlusEnd,1,M);
Ydiff = Y - Y'+eye(M);

WW = repmat(1./prod(Ydiff,2),1,M);
D = WW./(WW'.*Ydiff);

D(1:M1:M2) = 1-sum(D);
D = -D';                %Diff matrix
D = D(1:end-1,:);       %Augment for LGR points

Dd = [diag(diag(D(:,1:N1))), zeros(N1,1)]; % Diag diff matrix
Do = D - Dd; % off-diag diff matrix
