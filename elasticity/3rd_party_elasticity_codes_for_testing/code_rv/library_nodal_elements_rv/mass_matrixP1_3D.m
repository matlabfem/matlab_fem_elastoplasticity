function M=mass_matrixP1_3D(elements,volumes,coeffs)
%coeffs can be only P0 (elementwise constant) function 
%represented by a collumn vector with size(elements,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally
%Note: P1 coeffs needs a higher integration rule (not implemented yet)

if 1
    Xscalar=kron(ones(1,4),elements); Yscalar=kron(elements,ones(1,4)); 

    if (nargin<3)
        Zmassmatrix=kron(volumes,reshape((ones(4)+eye(4))/20,1,16)); 
    else
        if numel(coeffs)==size(elements,1) %P0 coefficients
            Zmassmatrix=kron(volumes.*coeffs,reshape((ones(4)+eye(4))/20,1,16)); 
        else %P1 coefficients
            M1=[6 2 2 2; 2 2 1 1; 2 1 2 1; 2 1 1 2]/120;
            M2=M1([4,1,2,3],[4,1,2,3]);
            M3=M2([4,1,2,3],[4,1,2,3]);
            M4=M3([4,1,2,3],[4,1,2,3]);

            Zmassmatrix=kron(volumes.*coeffs(elements(:,1)),reshape(M1,1,16)) ...
                       +kron(volumes.*coeffs(elements(:,2)),reshape(M2,1,16)) ...
                       +kron(volumes.*coeffs(elements(:,3)),reshape(M3,1,16)) ...
                       +kron(volumes.*coeffs(elements(:,4)),reshape(M4,1,16));                      
        end    
    end
    
    M=sparse(Xscalar,Yscalar,Zmassmatrix);
else
    NE=size(elements,1); %number of triangles    

    %particular part for a given element in a given dimension
    NLB=4; %number of local basic functions, it must be known!  
    IP=integration_point_transformation([-0.7236067977, -0.7236067977, -0.7236067977; ...
                                           0.1708203932, -0.7236067977, -0.7236067977; ...
                                          -0.7236067977,  0.1708203932, -0.7236067977; ...
                                          -0.7236067977, -0.7236067977,  0.1708203932]);
    weight =[1/4 1/4 1/4 1/4];
    phi=shapefun (IP,'P1');
    Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

    %copy this part for a creation of a new element
    M_local=zeros(NLB);
    for i=1:size(IP,2)
        M_local=M_local+weight(i)*phi(:,i)*phi(:,i)';
    end

    if (nargin<3)
        Z=astam(volumes,reshape(repmat(M_local,1,NE),NLB,NLB,NE));
    else
        Z=astam(volumes.*coeffs,reshape(repmat(M_local,1,NE),NLB,NLB,NE));
    end
    X=permute(Y,[2 1 3]);  
    M=sparse(X(:),Y(:),Z(:));
end