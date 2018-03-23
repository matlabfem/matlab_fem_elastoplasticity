function draw_vertices(coordinates, vertices)
if (nargin==1)
    vertices=1:size(coordinates,1);
end
if (size(vertices,2)>size(vertices,1))
    vertices=vertices';
end

plot(coordinates(vertices,1),coordinates(vertices,2),'o');
end
