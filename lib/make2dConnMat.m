function M=make2dConnMat(nX,addrow)

nMat=nX^2;
if addrow
    nMat=nMat+1;
end
whichDiags=[-nX -1 0 1 nX];

offDiags=-ones(1,nMat);
offDiags(mod(1:nMat,nX)==0)=0;

diagVals=[-ones(nMat,1) offDiags' 4*ones(nMat,1) offDiags', -ones(nMat,1)];
M=spdiags(diagVals,whichDiags,nMat,nMat);
