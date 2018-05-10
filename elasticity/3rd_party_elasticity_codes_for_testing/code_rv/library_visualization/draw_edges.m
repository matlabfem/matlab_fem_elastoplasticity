function draw_edges(edges,coordinates,format,Zmax,which)

if (nargin<5)
    which=1:size(edges,1);
end
if (size(which,2)>size(which,1))
    which=which';
end

if (nargin<4)
    Zmax=0;
end

X=[coordinates(edges(which,1),1) coordinates(edges(which,2),1)]';
Y=[coordinates(edges(which,1),2) coordinates(edges(which,2),2)]';
if (size(coordinates,2)==2)   
    plot3(X,Y,Zmax*ones(size(Y)),format,'LineWidth',0.5)
else
    Z=[coordinates(edges(which,1),3) coordinates(edges(which,2),3)]';
    plot3(X,Y,Z,format,'LineWidth',0.5)
end
