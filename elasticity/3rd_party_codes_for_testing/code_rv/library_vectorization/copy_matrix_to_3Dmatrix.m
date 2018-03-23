function Z = copy_matrix_to_3Dmatrix(Z_local,NE)
Z=zeros(size(Z_local,1),size(Z_local,2),NE);
for i=1:size(Z_local,1)
    for j=1:size(Z_local,2)
        Z(i,j,:)=Z_local(i,j);
    end
end

end

