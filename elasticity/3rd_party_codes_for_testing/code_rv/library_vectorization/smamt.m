function amb = smamt (smx,ama)
% ama: ama(1:ny,1:nx,1:nz)
% smx: smx(1:nk,1:nx)
% amb: amb(1:nk,1:ny,1:nz)

[ny,nx,nz] = size(ama);
[nk,nx]    = size(smx);

amb     = zeros(nk,ny,nz);
for row = 1:nk

    amb(row,:,:) = svamt(smx(row,:),ama);
    
end

return