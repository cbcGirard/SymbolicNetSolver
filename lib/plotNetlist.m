function plotNetlist(nodes,rTerms,rSize,rVals,rColors,darkMode)
if nargin<5
    rColors=repmat('k',size(rTerms,1),1);
    darkMode=false;
end
for ii=1:length(rVals)
    cent(ii,:)=mean(nodes(rTerms(ii,:),:));
end
for ii=1:length(rVals)
    theNodes=nodes(rTerms(ii,:),:);
    collision=false;
    for jj=1:ii-1
        distTo=norm(cent(ii,:)-cent(jj,:));
        collision=collision||(distTo<rSize(1));
        
    end
    plotResistor(collision,theNodes,rSize,rColors(ii,:),rVals(ii),darkMode)
    axis equal
    axis off
end


function plotResistor(collision,nodeLocs,rSize,colorspec,valname,darkMode)

if darkMode
    lineColor=[0.68627451, 0.811764706, 1];
    fillColor=[0, 0.098039216, 0.137254902];
else
    lineColor='k';
    fillColor='w';
end

wR=rSize(1);
lR=rSize(2);
vec=nodeLocs(1,:)-nodeLocs(2,:);
ang=atan(vec(1)/vec(2));

hcorners=[wR wR -wR -wR; lR -lR -lR lR]';
rotmat=[cos(ang) sin(ang); -sin(ang) cos(ang)];
corners=zeros(4,2);

center=mean(nodeLocs);
if collision
    center=center+(rotmat*[-3*wR; 0])';
end

for ii=1:4
    corners(ii,:) = center+(rotmat*hcorners(ii,:)')';
end

if collision
    nodeLocs=[nodeLocs(1,:); center; nodeLocs(2,:)];
end

plot(nodeLocs(:,1),nodeLocs(:,2),'Color',lineColor)

fill(corners(:,1),corners(:,2),NaN,'EdgeColor',colorspec,'FaceColor',fillColor)

if ~isempty(valname)
    if isa(valname,'sym')
        valname=char(valname);
    end
    text(center(1), center(2),valname,...
        'HorizontalAlignment','center','Color',colorspec)
end
