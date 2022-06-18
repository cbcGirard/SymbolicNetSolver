%% Setup

% create Tex document containing transform steps
fileOutput=true;
% resultFolder=pwd;
resultFolder='Y:\Diagrams, etc\Rmesh';

darkMode=true;

onDiag=false;
altConfig=false;

if fileOutput
    folder=[resultFolder '/output/edge'];
    if exist(folder,'dir')==0
        mkdir(folder);
    end
else
    folder=false;
end

if onDiag
    testNodes=[1 9];
else
    testNodes=[1 3];
end

%% Initialize

if fileOutput
    texfile=fopen([folder '/results.tex'],'w');
else
    texfile=1;
end
nX=2;
syms R_E R_C real;
[~,~,~,Rx,Ry]=make2dMixedConn(nX,[1,0],[R_E R_C/R_E]);
vLoc=[1 0];

if altConfig
    R=[ones(4,1)*R_E; ones(4,1)*R_C];
    nodes=sort([1 3; 3 9; 9 7; 7 1; 1 5; 7 5; 3 5; 9 5],2);
else
    [R,nodes]=toNetlist(nX,Rx,Ry);
end
xx=linspace(0,1,nX+1);
[XX YY]=meshgrid(xx,xx');
nodeLocs=[XX(:) YY(:)];

%% Along edge
fprintf(texfile,'\\subsection{Refined mesh, along edge}\\label{refine10}\n');
[Req, rvals]=simplifyRNet_v2(nodeLocs,R,nodes,[1 3],folder,texfile,darkMode);
% for ii=1:length(rvals)
%     if folder
%         fprintf(texfile,'$$ R_%i = %s $$\n',ii,latex(rvals(ii)));
%     else
%         fprintf('R_%i = %s\n',ii,rvals(ii));
%     end
% end
if folder
    fprintf(texfile,'\n\n $$ R_{eq} = %s $$ \n\n', latex(Req));
else
    fprintf('\n\n Req = %s\n\n',Req)
end

fclose(texfile)
%% Across corners
fprintf(texfile,'\\subsection{Refined mesh,across corners}\\label{refine11}\n');

[Req,rvals]=simplifyRNet_v2(nodeLocs,R,nodes,[1 9],folder,texfile);

if folder
    fprintf(texfile,'\n\n $$ R_{eq} = %s $$ \n\n', latex(Req));
else
    fprintf('\n\n Req = %s\n\n',Req)
end

%% Interpolate from nodes only
R=[ones(4,1)*R_E; ones(4,1)*R_C];
nodes=sort([1 3; 3 9; 9 7; 7 1; 1 5; 7 5; 3 5; 9 5],2);

%% Along edge
fprintf(texfile,'\\subsection{Interpolated mesh, along edge}\\label{interp10}\n');

[Req, rvals]=simplifyRNet_v2(nodeLocs,R,nodes,[1 3],folder,texfile);

if folder
    fprintf(texfile,'\n\n $$ R_{eq} = %s $$ \n\n', latex(Req));
else
    fprintf('\n\n Req = %s\n\n',Req)
end

%% Across corners
fprintf(texfile,'\\subsection{Interpolated mesh, across corners}\\label{interp11}\n');

[Req, rvals]=simplifyRNet_v2(nodeLocs,R,nodes,[1 9],folder,texfile);

if folder
    fprintf(texfile,'\n\n $$ R_{eq} = %s $$ \n\n', latex(Req));
else
    fprintf('\n\n Req = %s\n\n',Req)
end
