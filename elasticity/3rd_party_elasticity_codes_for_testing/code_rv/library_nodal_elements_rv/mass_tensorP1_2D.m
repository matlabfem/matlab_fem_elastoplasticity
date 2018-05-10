function Z=mass_tensorP1_2D(elements,areas,coeffs)
%coeffs can be only P0 (elementwise constant) function 
%represented by a collumn vector with size(elements,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally
%Note: P1 coeffs needs a higher integration rule (not implemented yet)

NE=size(elements,1); %number of elements
Z_local=(ones(3)+eye(3))/12; 

Z=copy_matrix_to_3Dmatrix(Z_local,NE);

if (nargin<3)
    Z=astam(areas',Z);
else
    if numel(coeffs)==size(elements,1) %P0 coefficients
        Z=astam((areas.*coeffs)',Z);
    end      
end