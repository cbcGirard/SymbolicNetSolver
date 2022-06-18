function setupMixed(gen,vLoc)

nX=2^gen;

iX=zeros(nX,nX+1);
iY=iX';
V=zeros(nX+1);

[M, v,edge]=make2dMixedConn(nX,vLoc);

rX=[edge(:,1); edge(end-nX+1:end,4)];
rY=[edge(:,2); edge(nX:nX:end,3)];

