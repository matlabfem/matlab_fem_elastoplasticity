function show_edges(edge2nodes,coordinates,which)
%function: show_edges(edge,coordinates,whichelements)
%requires: show_points
%displays edges numbers at their middle point coordinates  
%if which are not provided, all are displayed

if (nargin==2)
    %all elements will be displayed
    which=1:size(edge2nodes,1);
elseif (size(which,2)>size(which,1))
    which=which';
end

middle_points=(coordinates(edge2nodes(:,1),:)+coordinates(edge2nodes(:,2),:))/2;
if (size(coordinates,2)==3)
    %3D case
    show_points_blue(middle_points(which,1), middle_points(which,2),middle_points(which,3),which);
else
    %2D case   
    show_points_blue(middle_points(which,1), middle_points(which,2),0*middle_points(which,1),which);
end
