function [R iVec]=rMeshSym(nGen,atDiag)
syms R b M rational;

nX=2^(nGen)+1;
nLoops=4^nGen+1;
nR=(nX^2-nX);

%ind=@(x,y) x+y*(nX-1)+1;
ndxLoop=@(x,y) ind(x,y,nX-2,nX-2);


%I=zeros(
b=sym(zeros(nLoops,1));
b(end)=sym(1);


M=sym(4*eye(nLoops));
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
        M(curLoop,neighbors)=sym(-1);
        
        isedge=[jj==0 ii==(nX-2)];
        if any(isedge)
            if atDiag
                k=sym(-1);
                if all(isedge)
                    k=sym(-2);
                end
                
            else
                k=sym(-1);
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

iVec=linsolve(M,b);
R=1/iVec(end);

function ndx=ind(x,y,maxX,maxY)
isValid=(x>=0)&(x<=maxX)&(y>=0)&(y<=maxY);
ndx= x+(y)*(maxX+1)+1;
%ndx(~isValid)=[];
ndx=ndx(isValid);
ndx=unique(ndx);
