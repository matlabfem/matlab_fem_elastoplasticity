function show_constant_scalar(constantValue,nodes2coord,elems2nodes,nodalDisplacement)
    
if (nargin==4)     
    if norm(nodalDisplacement)==0
       factor=1;
    else
       factor =10^(-round(log10(max(max(nodalDisplacement)))));
    end   
else
    nodalDisplacement=zeros(size(nodes2coord));
end

X=nodes2coord(:,1)+nodalDisplacement(:,1);
Y=nodes2coord(:,2)+nodalDisplacement(:,2);
fill3(X(elems2nodes)',Y(elems2nodes)',kron(ones(1,size(elems2nodes,2)),constantValue)',kron(ones(1,size(elems2nodes,2)),constantValue)','FaceColor','interp','LineStyle','none');



end
