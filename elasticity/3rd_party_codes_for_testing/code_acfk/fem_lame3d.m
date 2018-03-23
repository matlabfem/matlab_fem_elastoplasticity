%FEM_LAME3D  three-dimensional finite element method for the Lame-Problem
%
%    FEM_LAME3D solves the Lame Problem
%      (lambda+mu)(grad div u) +mu div grad u +f = 0  in Omega
%                                            M u = w  on the Dirichlet boundary
%           (lambda tr eps(u) Id +2 mu eps(u)) n = g  on the Neumann boundary
%    on a geometry described by tetraeder and presents the solution graphically.
%
%    The parameters E and NU, defined in the first line, are problem
%    dependend and have to be chosen by the user.
%
%    Remark: This program is a supplement to the paper
%    "Matlab-Implementation of the Finite Element Method in Elasticity" by  
%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose. The reader should 
%    consult that paper for more information.   
%
%
%    M-files you need to run FEM_LAME3D
%       <stima.m>, <f.m>, <u_d.m>, <show.m>, and <g.m> (optional)
%
%    Data-files you need to run FEM_LAME3D
%       <coordinates.dat>, <elements.dat>, 
%       <dirichlet.dat>, <neumann.dat> (optional)

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <fem_lame3d.m> in $(HOME)/acfk/fem_lame3d/angle/
%    This program and corresponding data-files give Fig. 8 in
%    "Matlab-Implementation of the Finite Element Method in Elasticity"

% Initialisation
E=100000;nu=0.3;mu=E/(2*(1+nu));lambda=E*nu/((1+nu)*(1-2*nu));
load coordinates.dat;
load elements.dat;
eval('load neumann.dat;','neumann = [];');
load dirichlet.dat;

fprintf('Technique of Alberty, Carstensen, Funken, Klose: ')
tic    
%Assembly (original from paper)    
A = sparse(3*size(coordinates,1),3*size(coordinates,1));
for j = 1:size(elements,1)
    I = 3*elements(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
    A(I,I) = A(I,I) +stima(coordinates(elements(j,:),:),lambda,mu);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% added by Jan Valdman 
fprintf('Technique of Rahman and Valdman: ')
tic %Assembly of Rahman, Valdman for faster performance
A2=stiffness_matrixP1_3D_elasticity(elements,coordinates,lambda,mu);  
toc

fprintf('Technique of Cermak, Sysala and Valdman: ') 
tic %Assembly of Cermak, Sysala, Valdman for even faster performance   
shear = E/(2*(1+nu)) ;        % shear modulus
bulk = E/(3*(1-2*nu)) ;       % bulk modulus
[Xi, WF] = quadrature_volume_3D('P1');
[HatP,DHatP1,DHatP2,DHatP3] = local_basis_volume_3D('P1', Xi);
A3=stiffness_matrix_3D(elements',coordinates',shear,bulk,DHatP1,DHatP2,DHatP3,WF); 
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

b = zeros(3*size(coordinates,1),1);
%Volumeforces
for j = 1:size(elements,1)
  I = 3*elements(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0]; 
  fs = f(sum(coordinates(elements(j,:),:))/4)';
  b(I) = b(I) +det([1,1,1,1;coordinates(elements(j,:),:)'])*[fs;fs;fs;fs]/24;
end

%Neumann conditions
if ~isempty(neumann)
  for j = 1:size(neumann,1)
    n = cross( coordinates(neumann(j,2),:)-coordinates(neumann(j,1),:), ...
	coordinates(neumann(j,3),:)-coordinates(neumann(j,1),:));
    I = 3*neumann(j,[1,1,1,2,2,2,3,3,3])-[2,1,0,2,1,0,2,1,0];
    gm = g(sum(coordinates(neumann(j,:),:))/3,n/norm(n))';
    b(I) = b(I) +norm(n)*[gm;gm;gm]/6;   
  end
end

%Dirichlet conditions
dirichletnodes = unique(dirichlet);
[W,M] = u_d(coordinates(dirichletnodes,:));
B = sparse(size(W,1),3*size(coordinates,1));
for k = 0:2
  for l = 0:2
    B(1+l:3:size(M,1),3*dirichletnodes-2+k) = diag(M(1+l:3:size(M,1),1+k));
  end
end
mask=find(sum(abs(B)'));
A = [A,B(mask,:)';B(mask,:),sparse(length(mask),length(mask))];
b = [b;W(mask,:)];

%Calculation of the solution
x = A\b;
u = x(1:3*size(coordinates,1));

%Graphical representation
show(elements,dirichlet,neumann,coordinates,u,lambda,mu); 
