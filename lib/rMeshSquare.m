function [Req,V,Ix, Iy,M]=rMeshSquare(nGen,atDiag)

nX=2^(nGen)+1;
nLoops=4^nGen+1;
nR=(nX^2-nX);

%ind=@(x,y) x+y*(nX-1)+1;
ndxLoop=@(x,y) ind(x,y,nX-2,nX-2);

V=zeros(nX);
V(end,end)=1;
%I=zeros(
b=zeros(nLoops,1);
b(end)=1;


M=4*speye(nLoops);
if atDiag
M(end,end)=2^(nGen+1);
else
    M(end,end)=2^nGen;
end

for jj=0:nX-2
    for ii=0:nX-2
        curLoop=ndxLoop(ii,jj);
        
        %    if (ii-1)>=0
        %      M(curLoop,ind(ii-1,jj))=-1;
        %    end
        %   if (jj-1)>=0
        %      M(curLoop,ind(ii,jj-1))=-1;
        %    end
        %    if (ii+1)<=(nX-2)
        %      M(curLoop,ind(ii+1,jj))=-1;
        %    end
        %  if (jj+1)<=(nX-2)
        %      M(curLoop,ind(ii,jj+1))=-1;
        %    end
        nearX=ii+[0 1 0 -1];
        nearY=jj+[-1 0 1 0];
        neighbors=ndxLoop(nearX,nearY);
        M(curLoop,neighbors)=-1;
        
        isedge=[jj==0 ii==(nX-2)];
        if any(isedge)
            if atDiag
                k=-1;
                if all(isedge)
                    k=-2;
                end
                
            else
                k=-1;
                if jj~=0
                    break;
                end
            end
            M(end,ndxLoop(ii,jj))=k;
            M(ndxLoop(ii,jj),end)=k;
        end
        
        %    printf('ii=%d, jj=%d \n',ii,jj);
        %    printf('Loop %d: ',curLoop);
        %    printf('%d \t',M(curLoop,:));
        %    printf('\n');
        
    end
end

%full(M)

iVec=M\b;
Req=1/iVec(end);
%printf('%.2f\n',iVec);

iL=reshape(iVec(1:end-1),nX-1,nX-1);

if nGen==0
   xpadVal= 1;
   ypadVal=1;
else
    xpadVal=-M(end,1);
    ypadVal=-M(end,end-1);
end
xpad=iVec(end)*xpadVal*ones(nX+1,1);
ypad=iVec(end)*ypadVal*ones(1,nX-1);


iL=[ zeros(1,nX-1); iL; ypad];
iL=[ xpad iL zeros(nX+1,1)];


Iy=iL(2:end,2:end-1)-iL(1:end-1,2:end-1);
Ix=iL(2:end-1,2:end)-iL(2:end-1,1:end-1);

V=cumsum([cumsum([0 Iy(1,:)]);-Ix]);




function ndx=ind(x,y,maxX,maxY)
isValid=(x>=0)&(x<=maxX)&(y>=0)&(y<=maxY);
ndx= x+(y)*(maxX+1)+1;
%ndx(~isValid)=[];
ndx=ndx(isValid);
ndx=unique(ndx);
