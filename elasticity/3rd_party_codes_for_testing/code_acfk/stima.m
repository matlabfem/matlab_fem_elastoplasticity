function stima = stima(vertices,lambda,mu)
%STIMA   Computes element stiffness matrix for tetraeder.
%   M = STIMA(X,LAMBDA,MU) computes element stiffness matrix for
%   tetraeder. The coordinates of the vertices are stored in X. LAMBDA
%   and MU are the Lame constants.
%
%   This routine should not be modified.
%
%
%   See also FEM_LAME3D.

%    J. Alberty, C. Carstensen and S. A. Funken  07-03-00
%    File <stima.m> in $(HOME)/acfk/fem_lame3d/angle/

PhiGrad = [1,1,1,1;vertices']\[zeros(1,3);eye(3)];
R = zeros(6,12);
R([1,4,5],1:3:10) = PhiGrad';
R([4,2,6],2:3:11) = PhiGrad';
R([5,6,3],3:3:12) = PhiGrad';
C([1:3],[1:3]) = lambda*ones(3,3)+2*mu*eye(3);
C([4:6],[4:6]) = mu*eye(3);
stima = det([1,1,1,1;vertices'])/6*R'*C*R;
