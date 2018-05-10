function show_nodal_scalar(nodalValue,elements,coordinates,nodalDisplacement)
if nargin==3
    h=show_nodal_scalar_frame(nodalValue,elements,coordinates);
    set(h,'edgecolor','none')
elseif nargin==4
    h=show_nodal_scalar_frame(nodalValue,elements,coordinates,nodalDisplacement);
    set(h,'edgecolor','none')
end
    
