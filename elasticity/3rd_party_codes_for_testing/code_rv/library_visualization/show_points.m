function show_points(pointsX,pointsY,pointsZ,labels)
if (size(labels,2)>size(labels,1))
    labels=labels';
end
text(pointsX, pointsY, pointsZ,int2str(labels), 'fontsize', 14, 'fontweight', 'bold');

%'BackgroundColor',[0 0 1] - blue color