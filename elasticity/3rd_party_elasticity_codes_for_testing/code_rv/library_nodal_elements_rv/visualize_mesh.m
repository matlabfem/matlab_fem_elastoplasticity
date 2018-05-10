if (level<2)
    subplot(1,2,level+1); show_mesh(elements,coordinates); axis image; title(strcat('mesh, level=',num2str(level)));
    
%     if demo>0   
%         dn=unique(dirichlet); nn=unique(neumann);
%         hold on; show_nodes(coordinates,dn,nn); hold off    
%     end
end