function plotRNet(nX,vLoc)

xx=linspace(0,1,nX+1);
centeredR=@(cX,cY,dX,dY)[cX*ones(1,4); cY*ones(1,4)] +0.5*[-dX dX dX -dX; -dY -dY dY dY];
dx=diff(xx(1:2));
rL=2*dx/3;
rW=dx/4;

for ii=1:length(xx)
    plot(ones(1,2)*xx(ii),[0 1],'k',[0 1],ones(1,2)*xx(ii),'k');
end

for jj=1:length(xx)
    for ii=1:length(xx)
        if ii<length(xx)
            %horizontal
            cpoint=[mean(xx(ii:ii+1)) xx(jj)];
            corners=centeredR(cpoint(1),cpoint(2),rL,rW);
            fill(corners(1,:),corners(2,:),NaN,'EdgeColor','k','FaceColor','w');
        end
        
        %vertical
        if jj<length(xx)
            cpoint=[xx(ii) mean(xx(jj:jj+1))];
            corners=centeredR(cpoint(1),cpoint(2),rW,rL);
            fill(corners(1,:),corners(2,:),...
                NaN,'EdgeColor','k','FaceColor','w');
        end
    end
end

if nargin==2
    lstem=0.1;
    plot([0 0],[0 -lstem],'k')
    fill([-lstem lstem 0],[-lstem -lstem, -2*lstem],...
        NaN,'EdgeColor','k','FaceColor','w');
    
    if vLoc(2)==0
        plot(ones(1,2)*vLoc(1),...
            ones(1,2)*vLoc(2)-[0 lstem],'k')
        plot(ones(1,2)*vLoc(1)+[-lstem lstem],...
            -ones(1,2)*lstem,'k')
        text(vLoc(1),vLoc(2)-1.5*lstem,'\textbf{V+}',...
            'HorizontalAlignment','Center',...
            'Interpreter','latex')
    else
        plot(ones(1,2)*vLoc(1)+[0 lstem],...
            ones(1,2)*vLoc(2),'k')
        plot(ones(1,2)*(vLoc(1)+lstem),...
            ones(1,2)*vLoc(2)+[-lstem lstem],'k')
        text(vLoc(1)+lstem,vLoc(2),'\textbf{V+}',...
            'Interpreter','latex')
        
    end
    
end