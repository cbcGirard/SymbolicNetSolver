function [Mix,Miy]=make2dReconMat(nX,vLoc)

nRows=nX^2;

ndxVec=0:nRows-1;
isEdge=[ndxVec'<nX, mod(ndxVec',nX)==0,...
    mod(ndxVec',nX)==nX-1, ndxVec'>=nX*(nX-1)];

% Rxe=[isEdge(:,1); isEdge(end-nX+1:end,4)];
Rye=[isEdge(:,2); isEdge(nX:nX:end,3)];
nR=length(Rye);


Miy=spdiags([circshift(~Rye,-1),-ones(size(Rye))],...
    [-1 0],nR,nRows+1);
Miy(nRows+1,nRows+1)=0;
yext=0:nX-1<vLoc(2)*nX;
for ii=1:nX
    Miy(nX^2+ii,nX*ii)=1;
    Miy(nX^2+ii,end)=-yext(ii);
end

Mix=spdiags([ones(nR,1) -ones(nR,1)],[0 -nX],...
    nR,nRows+1);
Mix(nRows+1,nRows+1)=0;
Mix(0:nX-1<vLoc(1)*nX,end)=-1;