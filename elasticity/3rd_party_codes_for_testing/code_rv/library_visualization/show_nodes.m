function show_nodes(coordinates,dn,nn)

if (nargin==1)
    dn=1:size(coordinates,1);
    nn=[];
end


if size(coordinates,2)==2
    scatter(coordinates(dn,1),coordinates(dn,2),30,'o','MarkerFaceColor', 'b')
    hold on
    scatter(coordinates(nn,1),coordinates(nn,2),50,'ro')
    hold off
    xlabel('x'); ylabel('y');
elseif size(coordinates,2)==3
    %plot_points([coordinates(nn,1),coordinates(nn,2),coordinates(nn,3)],'rx'); 
    scatter3(coordinates(dn,1),coordinates(dn,2),coordinates(dn,3),30,'o', 'MarkerFaceColor', 'b')
    hold on
    scatter3(coordinates(nn,1),coordinates(nn,2),coordinates(nn,3),50, 'ro')
    hold off
    xlabel('x'); ylabel('y'); zlabel('z');
end
end

