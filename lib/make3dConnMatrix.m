function [M mmi]=make3dConnMatrix(nX,addrow)

nnormals=3;
nplanes=nX+1;
nloops=nX^2;

Nloops=nnormals*nplanes*nloops;
coords=zeros(Nloops,4);
for nn=0:nnormals-1
    maxn=(nX-1)*ones(1,3);
    maxn(nn+1)=nX;
    for kk=0:maxn(1)
        for jj=0:maxn(2)
            for ii=0:maxn(3)
%                 nearX
                priorLoops=1+nn*(nplanes*nloops);
                mth=ii+nX*jj+(nX^2)*kk;
                coords(mth,:)=[nn ii jj kk];
%                 disp(mth);
                mmi(mth)=mth;
            end
        end
    end
end
M=coords;