function [M,v,isEdge,Rx,Ry,Mix,Miy]=make2dMixedConn(nX,vLoc,ekVal)


nRows=nX^2;
ndxVec=0:nRows-1;
isEdge=[ndxVec'<nX, mod(ndxVec',nX)==0,...
    mod(ndxVec',nX)==nX-1, ndxVec'>=nX*(nX-1)];
neighbors=[-nX -1 1 nX];

edges=spalloc(nRows+1,nRows+1,4*nRows);
centers=edges;

for ii=0:length(ndxVec)-1
    nxt=neighbors(~isEdge(ii+1,:))+ii+1;
    centers(ii+1,nxt)=-1;
    nEdges=sum(isEdge(ii+1,:));
    centers(ii+1,ii+1)=4-nEdges;
    edges(ii+1,ii+1)=nEdges;
    
end

extConn=make2dEdgeConn(nX,vLoc);
edges(end,:)=extConn;
edges(:,end)=extConn';

if nargin==3
    E=ekVal(1);
    k=ekVal(2);
    me=E*edges;
    mc=E*k*centers;
else
syms E k;
me=sym(E*edges);
mc=sym(E*k*centers);
end
M=mc+me;

v=zeros(nRows+1,1);
v(end)=1;

Rxe=[isEdge(:,1); isEdge(end-nX+1:end,4)];
Rye=[isEdge(:,2); isEdge(nX:nX:end,3)];

Rx=E*Rxe+E*k*not(Rxe);
Ry=E*Rye +E*k*not(Rye);
nR=length(Rx);

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