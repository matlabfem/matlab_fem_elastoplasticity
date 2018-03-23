function dem = amdet (ama)
% ama: ama(1:nx,1:nx,1:nz)
% amb: amb(1:nx,1:nx,1:nz)

[nx,nx,nz] = size(ama);

if nx > 3, 
   
   msg = 'Array operation for invertig matrices of dimension more ';
   msg = [msg 'than 3 (three) is not available, so the performance may be bad.'];
   error(msg);
   
end

if nx==1
    
   % Matrix determinant.
   dem = squeeze(ama);
   
elseif nx==2

   % Matrix elements.
   a = squeeze(ama(1,1,:)); b = squeeze(ama(1,2,:));
   c = squeeze(ama(2,1,:)); d = squeeze(ama(2,2,:));

   % Matrix determinant.
   dem = a.*d - b.*c;
   
else 
    % Matrix elements.
   a = squeeze(ama(1,1,:)); b = squeeze(ama(1,2,:)); c = squeeze(ama(1,3,:)); 
   d = squeeze(ama(2,1,:)); e = squeeze(ama(2,2,:)); f = squeeze(ama(2,3,:));  
   g = squeeze(ama(3,1,:)); h = squeeze(ama(3,2,:)); i = squeeze(ama(3,3,:));   
    
   % Matrix determinant.
   dem = a.*(e.*i - f.*h)-b.*(d.*i-f.*g)+c.*(d.*h-e.*g);
end

dem = dem(:)';

return