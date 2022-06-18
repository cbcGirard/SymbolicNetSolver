function [rvals, rAll]=simplifyRNet_v2(netNodes,Rvals,Rterms,fixedNodes,texfolder,texfile,darkMode)
% Master program to reduce symbolic netlist to equivalent resistance
% between nodes

if nargin<7
    darkMode=false;
end
global rterms rvals rAll rNdx nodelocs
rterms=Rterms;
rvals=Rvals;
rAll=unique(rvals);
q=double(repmat(rvals,1,length(rAll))==repmat(rAll',length(rvals),1));
rNdx=q*(1:length(rAll))';

nodelocs=netNodes;

netBnds=[[min(netNodes(:,1)), min(netNodes(:,2))];
        [max(netNodes(:,1)), max(netNodes(:,2))]];
netCtr=0.5*(netBnds(1,:)+netBnds(2,:));
netSpan=0.5*(netBnds(2,:)-netBnds(1,:));
bndPts=[netCtr-1.5*netSpan, netCtr+1.5*netSpan];

figure;
 hold all;
 axis manual;
 axis(bndPts([1,3,2,4]));
 set(gca)
 if darkMode
     scope=importdata("lib\cmapScope.mat");
     cmap=scope;

 else
cmap=colormap(lines(length(unique(rNdx))));
 end

rnames=strsplit(sprintf('R_{%i} ',rNdx));
plotNetlist(nodelocs,rterms,[0.05 0.15],rnames(1:end-1),cmap(rNdx,:),darkMode)
if texfolder
    %     [~,fname]=fileparts(tempname);
    %     fname=[isTex '/' fname];
    nfiles=length(dir(texfolder));
    fname=sprintf('fig%i',nfiles);
    print([texfolder '/' fname],'-depsc2');
    fprintf(texfile,'\\includegraphics[width=0.5\\textwidth]{%s}\n',fname);
else
    snapnow
end

while size(rterms,1)>length(fixedNodes)
    
    while any(countConn(rterms,fixedNodes)==2)
        %series resistors
        %         subplot 121
        %         plotNet(rvals,rterms,nodelocs)
        
        nconn=countConn(rterms,fixedNodes);
        elimNodes=find(nconn==2);
        elimR=find(any(rterms==elimNodes(1),2));
        newR=sum(rvals(elimR));
        if isa(newR,'sym')
            newR=simplify(newR);
        end
        newTerm=setxor(rterms(elimR(1),:),rterms(elimR(2),:));
        printupdate(elimR,newR,newTerm,'series',texfolder,texfile,darkMode)
    end
    if length(rvals)==1
        break
    end
    while size(isParallel(rterms),2)>1
        
        elimR=isParallel(rterms);
        elimR=elimR(1,:)';
        newR=parsum(rvals(elimR));
        newTerm=rterms(elimR(1),:);
        
        printupdate(elimR,newR,newTerm,'parallel',texfolder,texfile,darkMode)
    end
    if length(rvals)==1
        break
    end
    if any(countConn(rterms,fixedNodes)==3)
        
        nconn=countConn(rterms,fixedNodes);
        elimNodes=find(nconn==3);
        elimR=find(any(rterms==elimNodes(1),2));
        [newR, newTerm]=y2d(rvals(elimR),rterms(elimR,:));
        printupdate(elimR,newR,newTerm,'y2d',texfolder,texfile,darkMode)
    end
    
end

function printupdate(elimR,newR,newTerm,opType,isTex,texfile,darkMode)
global rNdx rterms nodelocs rvals rAll

switch opType
    case 'series'
        str1=sprintf('R_{%i} + R_{%i} = ',rNdx(elimR));
    case 'parallel'
        str1=sprintf('R_{%i} || R_{%i} = ',rNdx(elimR));
    case 'y2d'
        str1=sprintf('Y:[ R_{%i}, R_{%i}, R_{%i},] =\n',rNdx(elimR));
    case 'd2y'
        
end

updateNet(elimR,newR,newTerm);
switch opType
    case 'y2d'
        newNdx=rNdx(end-2:end);
        str2=sprintf('Delta: [ R_{%i}, R_{%i}, R_{%i},]',newNdx);
    case 'd2y'
        
    otherwise
        newNdx=rNdx(end);
        str2=sprintf('R_{%i}',newNdx);
end
if isTex
    fprintf(texfile,'$$ %s',str1);
    fprintf(texfile,'%s $$ \n',str2);
    for ii=1:length(newNdx)
        fprintf(texfile,'$$ R_{%i} = %s $$ \n',...
            newNdx(ii),latex(newR(ii)));
    end
else
    fprintf('%s \n',str1);
    fprintf('%s \n',str2);
    for ii=1:length(newNdx)
    fprintf('R_{%i} = %s \n',...
        newNdx(ii),newR(ii));
    end
end

if darkMode
    cmap=importdata('lib\cmapScope.mat');
else
    cmap=colormap(lines(length(rAll)));
end
rnames=strsplit(sprintf('R_{%i} ',rNdx));

cla;
plotNetlist(nodelocs,rterms,[0.05 0.15],rnames(1:end-1),...
    cmap(cycleIndices(rNdx,length(cmap)),:), darkMode)
axis off;

if isTex
    %     [~,fname]=fileparts(tempname);
    %     fname=[isTex '/' fname];
    nfiles=length(dir(isTex));
    fname=sprintf('fig%i',nfiles);
    fpath=[isTex '/' fname];
    print(fpath,'-vector','-depsc');
        if ispc
            print(fpath, '-vector','-dmeta');
        end
    print(fpath,'-dpng');
    %         pause(3);
    fprintf(texfile,'\\includegraphics[width=0.5\\textwidth]{%s}\n',fname);
else
    snapnow
end



function nuInd=cycleIndices(inds,nVals)
% prevent index overflow by cycling to beginning
nuInd=mod(inds-1,nVals)+1;

function updateNet(elimR,newR,newTerm)
global rterms rvals rAll rNdx;
rterms(elimR,:)=[];
rvals(elimR,:)=[];
rterms=[rterms;newTerm];
rvals=[rvals;newR];
for ii=1:length(newR)
    if any(newR(ii)==rAll)
        samendx=find(newR(ii)==rAll);
        rNdx=[rNdx; samendx];
    else
        rNdx=[rNdx; length(rAll)+1];
        rAll=[rAll;newR(ii)];
    end
end
rNdx(elimR)=[];

function nconn=countConn(rterms,fixedNodes)
nconn=zeros(max(rterms(:)),1);
for ii=1:length(nconn)
    if any(ii==fixedNodes)
    else
        nconn(ii)=sum(rterms(:)==ii);
    end
end


function [rY, termY]=d2y(rDelt,termDelt)
rd=sum(rDelt);
rY=sum(prod(nchoosek(rDelt,2),2))/rd;
pairs=nchoosek(1:length(rDelt),2);
for ii=1:length(pairs)
    termY=intersect(termDelt(pairs(ii,1),:),...
        termDelt(pairs(ii,2),:));
end
if isa(rY,'sym')
    rY=simplify(rY);
end

function [rDelt, termDelt]=y2d(rY,termY)
nR=length(rY);
if isa(rY,'sym')
    rDelt=sym(zeros(nR,1));
else
    rDelt=zeros(nR,1);
end
termDelt=zeros(nR,2);

rp=sum(prod(nchoosek(rY,2),2));
pairs=nchoosek(1:length(rY),2);
for ii=1:length(pairs)
    termDelt(ii,:)=setxor(termY(pairs(ii,1),:),...
        termY(pairs(ii,2),:));
    rDelt(ii)=rp/rY(setxor(1:nR,pairs(ii,:)));
end
if isa(rDelt,'sym')
    rDelt=simplify(rDelt);
end


function req=parsum(rPar)
req=1./(sum(1./rPar));
if isa(req,'sym')
    req=simplify(req);
end

function plotNet(rval,rterms,nodelocs)
% global rAll;
cla;
% figure
hold on
eachval=unique(rval);
cmap=colormap(lines(length(eachval)));
legob=zeros(length(eachval),1);
legName=cell(length(eachval),1);
% ptypes={'-','--','-.'};


%check for parallel resistors
sharesNodes=isParallel(rterms);

for ii=1:length(eachval)
    plotwhich=rval==eachval(ii);
    plotR=rval(plotwhich);
    plotT=rterms(plotwhich,:);
    nthPlace=find(plotwhich);
    for jj=1:length(plotR)
        %         ptype=ptypes{sharesNodes==jj};
        if any(nthPlace(jj)==sharesNodes(:,2:end))
            q=plot(nodelocs(plotT(jj,:),1),...
                nodelocs(plotT(jj,:),2),...
                '--','Color',cmap(ii,:),...
                'MarkerEdgeColor','k',...
                'Marker','.','LineWidth',4);
        else
            q=plot(nodelocs(plotT(jj,:),1),...
                nodelocs(plotT(jj,:),2),...
                'Color',cmap(ii,:),...
                'MarkerEdgeColor','k',...
                'Marker','.','LineWidth',2);
        end
        if jj==1
            legob(ii)=q;
            %             legName{ii}=char(eachval(ii));
            legName{ii}=sprintf('R_{%i}',eachval(ii));
        end
    end
end
legend(legob,legName,'Location','bestoutside')
xlim([-0.2 1.2])
ylim([-0.2 1.2])

axis equal

function parPairs=isParallel(rterm)
parPairs=[];
try
    for ii=1:length(rterm)
        whichPar=find(ismember(rterm,rterm(ii,:),'rows'))';
        if length(whichPar)>1
            %             break;
            parPairs=[parPairs; whichPar];
        end
    end
    parPairs=unique(parPairs,'rows');
catch
    parPairs=[];
end