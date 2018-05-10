function show_points_blue(pointsX,pointsY,pointsZ,labels)
if (size(labels,2)>size(labels,1))
    labels=labels';
end
text(pointsX, pointsY, pointsZ,int2str(labels),'color','b','fontsize', 10, 'fontweight', 'bold');