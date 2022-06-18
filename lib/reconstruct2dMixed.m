function [V, Ix, Iy]=reconstruct2dMixed(nX,I,Mix,Miy,Rx,Ry)
if nargin<6
    nR=size(Mix,1);
    Rx=ones(nR,1);
    Ry=Rx;
end

vIx=Mix*I;
vIy=Miy*I;

RRx=v2mx(Rx,nX);
RRy=v2my(Ry,nX);
Ix=v2mx(vIx,nX);
Iy=v2my(-vIy,nX);

Vx=Ix.*RRx;
Vy=Iy.*RRy;

Ve=cumsum([0 -Vy(1,:)]);
V=-cumsum([Ve;Vx]);

function M=v2mx(v,nX)
M=reshape(v,nX,nX+1);

function M=v2my(v,nX)
M=[reshape(v(1:nX^2),nX,nX); (v(nX^2+1:end)')];