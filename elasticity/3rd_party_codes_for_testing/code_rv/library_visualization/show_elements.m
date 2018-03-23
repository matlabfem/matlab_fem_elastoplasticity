function show_elements(elements,coordinates,which)
%function: show_elements(elements,coordinates,whichelements)
%requires: show_points
%displays element numbers at their middle point coordinates  
%if which are not provided, all are displayed

if (nargin==2)
    %all elements will be displayed
    which=1:size(elements,1);
elseif (size(which,2)>size(which,1))
    which=which';
end

if (size(coordinates,2)==3)
    %3D case
    middle_points=(coordinates(elements(:,1),:)+coordinates(elements(:,2),:)+coordinates(elements(:,3),:)+coordinates(elements(:,4),:))/4;
    show_points_red(middle_points(which,1), middle_points(which,2),middle_points(which,3),which);
else
    %2D case
    middle_points=evaluate_elements_average(elements,coordinates);
    show_points_red(middle_points(which,1), middle_points(which,2),0*middle_points(which,1),which);
end
