function circVec=make2dEdgeConn(nX,vLoc)

nMat=nX^2+1;

circVec=zeros(1,nMat);
circVec(1:nX*vLoc(1))=-1;

if vLoc(2)~=0
   circVec(nX:nX:(nMat-1)*vLoc(2))=-1;
   circVec(nX)=-2;
end
circVec(end)=-sum(circVec);
M=make2dConnMat(nX,1);
M(:,end)=circVec;
M(end,:)=circVec;

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
