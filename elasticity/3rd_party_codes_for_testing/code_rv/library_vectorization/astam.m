function amb = astam (asx,ama)
% ama: ama(1:nx,1:ny,1:nz)
% asx: avx(1,1:nz)
% amb: avb(1:nx,1:ny,1:nz)

[nx,ny,nz] = size(ama);

asx = reshape(asx,1,1,nz);
asx = asx(ones(nx,1),ones(ny,1),:);

amb = ama .* asx;

return