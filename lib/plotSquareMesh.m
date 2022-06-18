function plotSquareMesh(nX,V,Ix,Iy,plotI)

% [Req, V, Ix, Iy,M]=rMeshSquare(nGen,atDiag);

% syms R0 v0;
% % R0=sym(1/Req);
% sM=sym(M);
% sb=sym([zeros(size(M,1)-1,1); 1]);
% 
% sI=linsolve(sM,sb);
% R0=sym(1/sI(end));
% % fprintf(['Net current=' rats(v0/sI(end)) '\n']);
% fprintf(['Each resistor= ' char(R0) ' ohm \n']);

xx=linspace(0,1,nX+1);

% figure();
imagesc(xx,xx,V')
colormap jet
set(gca,'ydir','normal')
c=colorbar;
c.Label.String='Potential [V]';
xlabel('X [distance]');
ylabel('Y [distance]');
% title([sprintf('nX=%d, R_{eq}=',nX),rats(Req,25)]);



if plotI
    hold on
    plotRNet(xx);
    plotIvec(xx,Ix,Iy);
end

function plotIvec(xx,Ix,Iy)
dx=diff(xx(1:2));
cx=xx(2:end)-dx/2;
% cx=xx(2:end);

quiver(cx,xx,Ix',zeros(size(Ix')),0.5,'m');
quiver(xx,cx,zeros(size(Iy')),-Iy',0.5,'m');
