function amb = amtam (amx,ama)
% ama: ama(1:nx,1:ny,1:nz)
% amx: amx(1:nx,1:nk,1:nz)
% amb: amb(1:nk,1:ny,1:nz)

[nx,ny,nz] = size(ama);
[nx,nk,nz] = size(amx);

amb     = zeros(nk,ny,nz);
for row = 1:nk

    amb(row,:,:) = avtam(amx(:,row,:),ama);
    
end

return