function labelRnet(nX,Rx,Ry)
Ry2=[reshape(Ry(1:nX^2),nX,nX)' Ry(nX^2+1:end)];
Ry=Ry2(:);
xx=linspace(0,1,nX+1);
cx=diffcent(xx);
[XX CX]=meshgrid(xx,cx);

for ii=1:length(Rx)
    text(XX(ii),CX(ii),char(Ry(ii)),'HorizontalAlignment','center')
    text(CX(ii),XX(ii),char(Rx(ii)),'HorizontalAlignment','center')
end