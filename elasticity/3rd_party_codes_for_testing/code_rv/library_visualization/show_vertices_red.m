function show_vertices_red(coordinates, vertices)
if (nargin==1)
    vertices=1:size(coordinates,1);
end

if (size(vertices,2)>size(vertices,1))
    vertices=vertices';
end

show_points_red(coordinates(vertices,1),coordinates(vertices,2),0*coordinates(vertices,2),vertices);  %2D only modification
end
