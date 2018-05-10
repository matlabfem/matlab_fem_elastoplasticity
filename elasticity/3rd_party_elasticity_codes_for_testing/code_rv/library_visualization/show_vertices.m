function show_vertices(coordinates, vertices)
if (nargin==1)
    vertices=1:size(coordinates,1);  %all point will be displayed
end
if (size(vertices,2)>size(vertices,1))
    vertices=vertices';              %transposing the vector
end

if size(coordinates,2)==2
    show_points(coordinates(vertices,1),coordinates(vertices,2),0*coordinates(vertices,2),vertices);  %2D 
elseif size(coordinates,2)==3
    show_points(coordinates(vertices,1),coordinates(vertices,2),coordinates(vertices,3),vertices);  % 3D
end
