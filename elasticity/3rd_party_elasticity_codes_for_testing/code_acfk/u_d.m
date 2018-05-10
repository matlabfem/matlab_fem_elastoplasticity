function [W,M] = u_d(x)
%U_D   Data on the Dirichlet boundary
%   [W,M] = U_D(X) returns at N discrete points on the Dirichlet boundary
%   the direction for which the dispacement is given and the corresponding
%   values. If U is the displacement vector the Dirichlet boundary condition
%   is given by M*U = W. This input data has to be choosen by the user. X
%   has dimension N x 3, W has dimension 3*N x 1, and M has dimension 3*N
%   x 3.
%
%
%   See also FEM_LAME3D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <u_d.m> in $(HOME)/acfk/fem_lame3d/angle/

M = zeros(3*size(x,1),3);
W = zeros(3*size(x,1),1);
M(1:3:end,1) = 1;
M(2:3:end,2) = 1;
M(3:3:end,3) = 1;
W(3*find(x(:,2)<-50)) = 0.1;
