function M=mass_matrixP1_3D_elasticity(elements,areas)

%for elasticity in 2D
NE=size(elements,1); %number of elements

elements_elasticity=2*elements(:,[1 1 2 2 3 3])-kron(ones(NE,1),[1,0,1,0,1,0]);
NLB=6;
Y_3D=reshape(repmat(elements_elasticity,1,NLB)',NLB,NLB,NE);
X_3D=permute(Y_3D,[2 1 3]);

Zlocal=zeros(6);
Zlocal([1 3 5],[1 3 5])=(ones(3)+eye(3))/12;
Zlocal([2 4 6],[2 4 6])=(ones(3)+eye(3))/12;
Z_3D=astam(areas',repmat(Zlocal,[1 1 NE]));

M=sparse(X_3D(:),Y_3D(:),Z_3D(:));


%for laplace in 3D
NE=size(elements,1); %number of triangles    

%particular part for a given element in a given dimension
NLB=4; %number of local basic functions, it must be known!  
IP=integration_point_transformation([-0.7236067977, -0.7236067977, -0.7236067977; ...
                                       0.1708203932, -0.7236067977, -0.7236067977; ...
                                      -0.7236067977,  0.1708203932, -0.7236067977; ...
                                      -0.7236067977, -0.7236067977,  0.1708203932]);
weight =[1/4 1/4 1/4 1/4];
phi=shapefun (IP,'P1');
Y_3Dmatrix=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

%copy this part for a creation of a new element
M_local=zeros(NLB);
for i=1:size(IP,2)
    M_local=M_local+weight(i)*phi(:,i)*phi(:,i)';
end     
M_local_3Dmatrix=astam(areas,reshape(repmat(M_local,1,NE),NLB,NLB,NE));
X_3Dmatrix=permute(Y_3Dmatrix,[2 1 3]);  
M=sparse(X_3Dmatrix(:),Y_3Dmatrix(:),M_local_3Dmatrix(:));