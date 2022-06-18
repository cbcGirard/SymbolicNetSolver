function [R, nodes]=toNetlist(nX,Rx,Ry)
% Converts list of horizontal and vertical resistors into a netlist, 
% arranged in a nX by nX grid
nX=nX+1;
nnodes=1:nX^2;
R=[Rx; Ry];

xnodes=nnodes(mod(nnodes,nX)~=0);
sharednodes=xnodes(1:end-nX+1);
xedge=nnodes(end-nX+1:end);
yedge=nnodes(mod(nnodes,nX)==0);
nodes=[sharednodes, xedge(1:end-1), sharednodes, yedge(1:end-1);...
    sharednodes+1,xedge(2:end), sharednodes+nX,yedge(2:end)]';