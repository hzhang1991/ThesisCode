function CreateMatrixGlobal(morder, bctype, showplots)
%CREATEMATRIXGLOBAL    Create the Matrix for Poisson solves
%
% This function creates the matrix used for solving Poisson's equation.  It
% uses a continuous Galerkin formulation to define a potential, and then
% gradients automatically become discontinuous after differentiating the first
% formulation.
%
% Each application is expected to call this global function.
%
% When creating the mesh, the SubFactor parameter located in "input2D.data" 
% needs to match morder, the order of the method.
%
% Inputs:
%
% morder    = 1 (CG1 method)
%           = 2 (CG2 method -- needs to match SubFactor parameter)
%           = 3 (CG2 method -- needs to match SubFactor parameter)
%
% bctype    = 1 (default: Dirichlet BCs)
%           = 2 (double periodic BCs)
%
% showplots = 1 (default: show sparsity plots of matrices)
%             0 (don't show sparsity plots)
%
% Outputs:
%
% None, save for the following files that get written to 
% outputdir = './output':
%
%    Rnz.dat - number of nonzeros on each row
%    R.dat   - R matrix
%    Lnz.dat - number of nonzeros per row for L matrix
%    L.dat   - L matrix
%    P.dat   - P vector

outputdir = '.';
meshdir = ['../Unstructured_Mesh/mesh_output'];

if ( (morder~=1)&&(morder~=2)&&(morder~=3) )
    error([' Supported values of morder are 1, 2, and 3.  morder = ',num2str(morder)]);
end

fids = fopen([meshdir '/mesh_params.dat'],'r');
if fids==-1
    disp(' ');
    error(['File  ',meshdir,'/mesh_params.dat  not found.']);
end
NumElems      = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumPhysElems  = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumGhostElems = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumNodes      = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumPhysNodes  = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumBndNodes   = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
NumEdges      = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
fclose(fids);

fids = fopen([meshdir '/mesh_tnode.dat'],'r');
tmp = fscanf(fids,'%d %d %d',[3,inf]);
fclose(fids);
tnode = transpose(tmp);
clear tmp;

fids = fopen([meshdir '/mesh_node.dat'],'r');
tmp = fscanf(fids,'%e',[2,inf]);
fclose(fids);
node = transpose(tmp);
clear tmp;

fids = fopen([meshdir '/mesh_bnd_node.dat'],'r');
tmp = fscanf(fids,'%d',[1,inf]);
fclose(fids);
bnd_node = transpose(tmp);
clear tmp;

fids = fopen([meshdir '/mesh_area_prim.dat'],'r');
tmp = fscanf(fids,'%e',[1,inf]);
fclose(fids);
area = transpose(tmp);
clear tmp;

fids = fopen([meshdir '/mesh_jmat.dat'],'r');
Jmat = zeros(NumElems,2,2);
for i=1:NumElems
    tmp = fscanf(fids,'%e',[1,4]);
    Jmat(i,1,1) = tmp(1);
    Jmat(i,1,2) = tmp(2);
    Jmat(i,2,1) = tmp(3);
    Jmat(i,2,2) = tmp(4);
end
fclose(fids);
clear tmp;

xmin = min(node(1:NumPhysNodes,1));
xmax = max(node(1:NumPhysNodes,1));

ymin = min(node(1:NumPhysNodes,2));
ymax = max(node(1:NumPhysNodes,2));

if (morder>1)
    fids = fopen([meshdir '/submesh_params.dat'],'r');
    if fids==-1
        disp(' ');
        error(['File  ',meshdir,'/submesh_params.dat  not found.']);
    end
    SubFactor       = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    SubNumPhysElems = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    SubNumPhysNodes = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    SubNumBndNodes  = fscanf(fids,'%d',1);  fscanf(fids,'%s',1); fscanf(fids,'%s',1);
    fclose(fids);

    fids = fopen([meshdir '/submesh_tnode.dat'],'r');
    tmp = fscanf(fids,'%d %d %d',[3,inf]);
    fclose(fids);
    sub_tnode = transpose(tmp);
    clear tmp;
    
    fids = fopen([meshdir '/submesh_node.dat'],'r');
    tmp = fscanf(fids,'%e',[2,inf]);
    fclose(fids);
    sub_node = transpose(tmp);
    clear tmp;
    
    fids = fopen([meshdir '/submesh_bnd_node.dat'],'r');
    tmp = fscanf(fids,'%d',[1,inf]);
    fclose(fids);
    sub_bnd_node = transpose(tmp);
    clear tmp;

    if (morder==2)
        fids = fopen([meshdir '/submesh_node_subs.dat'],'r');
        tmp = fscanf(fids,'%d',[6,inf]);
        fclose(fids);
        node_subs = transpose(tmp);
        clear tmp;
    elseif (morder==3)
        fids = fopen([meshdir '/submesh_node_subs.dat'],'r');
        tmp = fscanf(fids,'%d',[10,inf]);
        fclose(fids);
        node_subs = transpose(tmp);
        clear tmp;
    end
end


% Construct stiffness matrix
if (bctype==1)
    
    if (morder==1)
        A = ConstructA_CG1(NumPhysNodes,...
                           NumPhysElems,...
                           NumBndNodes,...
                           node,...
                           tnode,...
                           area,...
                           bnd_node);
        SubNumBndNodes = NumBndNodes;
        sub_bnd_node = bnd_node;
        sub_node = node;
    elseif (morder==2)
        if (SubFactor~=2)
            error([' If morder=2, then need SubFactor=2,  SubFactor = ',num2str(SubFactor)]);
        end
        A = ConstructA_CG2(SubNumPhysNodes,...
                           NumPhysElems,...
                           SubNumBndNodes,...
                           sub_node,...
                           tnode,...
                           node_subs,...
                           area,...
                           sub_bnd_node,...
                           Jmat);
    elseif (morder==3)
        if (SubFactor~=3)
            error([' If morder=3, then need SubFactor=3,  SubFactor = ',num2str(SubFactor)]);
        end
        A = ConstructA_CG3(SubNumPhysNodes,...
                           NumPhysElems,...
                           SubNumBndNodes,...
                           sub_node,...
                           tnode,...
                           node_subs,...
                           area,...
                           sub_bnd_node,...
                           Jmat);
    end

elseif (bctype==2)
    
    if morder==1
        fids = fopen([meshdir '/mesh_node_map.dat'],'r');
        tmp = fscanf(fids,'%d',[1,inf]);
        fclose(fids);
        node_map = transpose(tmp);
        clear tmp;
    
        fids = fopen([meshdir '/mesh_rdnt_node.dat'],'r');
        NumRdntBndNodes = fscanf(fids,'%d',[1,1]);
        tmp = fscanf(fids,'%d',[2,inf]);
        fclose(fids);
        rdnt_bnd_node = transpose(tmp);
        clear tmp;
    
        SubNumRdntBndNodes = NumRdntBndNodes;
    else
        fids = fopen([meshdir '/submesh_node_map.dat'],'r');
        tmp = fscanf(fids,'%d',[1,inf]);
        fclose(fids);
        node_map = transpose(tmp);
        clear tmp;
        
        fids = fopen([meshdir '/submesh_rdnt_node.dat'],'r');
        SubNumRdntBndNodes = fscanf(fids,'%d',[1,1]);
        tmp = fscanf(fids,'%d',[2,inf]);
        fclose(fids);
        rdnt_bnd_node = transpose(tmp);
        clear tmp;
    end
    
    if (morder==1)
        A = ConstructAperiodic_CG1(NumPhysNodes,...
                                   NumPhysElems,...
                                   NumBndNodes,...
                                   NumRdntBndNodes,...
                                   node,...
                                   tnode,...
                                   area,...
                                   bnd_node,...
                                   node_map,...
                                   rdnt_bnd_node);
    elseif (morder==2)
        if (SubFactor~=2)
            error([' If morder=2, then need SubFactor=2,  SubFactor = ',num2str(SubFactor)]);
        end
        A = ConstructAperiodic_CG2(SubNumPhysNodes,...
                                   NumPhysElems,...
                                   SubNumBndNodes,...
                                   SubNumRdntBndNodes,...
                                   sub_node,...
                                   tnode,...
                                   node_subs,...
                                   area,...
                                   sub_bnd_node,...
                                   Jmat,...
                                   node_map,...
                                   rdnt_bnd_node);
    elseif (morder==3)
        if (SubFactor~=3)
            error([' If morder=3, then need SubFactor=3,  SubFactor = ',num2str(SubFactor)]);
        end
        A = ConstructAperiodic_CG3(SubNumPhysNodes,...
                                   NumPhysElems,...
                                   SubNumBndNodes,...
                                   SubNumRdntBndNodes,...
                                   sub_node,...
                                   tnode,...
                                   node_subs,...
                                   area,...
                                   sub_bnd_node,...
                                   Jmat,...
                                   node_map,...
                                   rdnt_bnd_node);
    end
    
end

P = amd(A);
Amod = A(P,P);
Rmod = chol(Amod);
Lmod = transpose(Rmod);

Amod_nz = nnz(Amod);
Rmod_nz = nnz(Rmod);
Lmod_nz = nnz(Lmod);

if morder==1
    ntmp = NumPhysNodes;
else
    ntmp = SubNumPhysNodes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output periodic boundary condition type
fid = fopen([outputdir '/boundary_condition.dat'],'w');
fprintf(fid,'%8i\n',bctype);
fclose(fid);

% Output subfactor
fid = fopen([outputdir '/subfactor.dat'],'w');
fprintf(fid,'%8i\n',morder);
fclose(fid);

% Output number of nonzeros on each row
Rmod_nz_per_row = zeros(ntmp,1);
for i=1:ntmp
    Rmod_nz_per_row(i) = nnz(Rmod(i,:));
end

% Output number of nonzeros per row for R matrix
fid = fopen([outputdir '/Rnz.dat'],'w');
for i=1:ntmp
    fprintf(fid,'%8i\n',Rmod_nz_per_row(i));
end
fclose(fid);

% Output R matrix
[index1,index2,Rvals]  = find(Rmod);
fid = fopen([outputdir '/R.dat'],'w');
for i=1:Rmod_nz    
    fprintf(fid,'%8i %8i %32.16e\n',index1(i),index2(i),Rvals(i));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output number of nonzers on each row
Lmod_nz_per_row = zeros(ntmp,1);
for i=1:ntmp
    Lmod_nz_per_row(i) = nnz(Lmod(i,:));
end

% Output number of nonzeors per row for L matrix
fid = fopen([outputdir '/Lnz.dat'],'w');
for i=1:ntmp
    fprintf(fid,'%8i\n',Lmod_nz_per_row(i));
end
fclose(fid);

% Output L matrix
[index1,index2,Lvals]  = find(Lmod);
fid = fopen([outputdir '/L.dat'],'w');
for i=1:Lmod_nz    
    fprintf(fid,'%8i %8i %32.16e\n',index1(i),index2(i),Lvals(i));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output P vector
fid = fopen([outputdir '/P.dat'],'w');
for i=1:ntmp
    fprintf(fid,'%8i\n',P(i));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (showplots==1)
    
    figure(1);
    clf;
    spy(A,5);
    set(gca,'fontsize',16);
    xlabel('');
    g1=text(0.05*ntmp,0.8*ntmp,['Non-zeros = ',num2str(Amod_nz)]);
    set(g1,'fontsize',16);
    set(g1,'FontWeight','bold');
    title(['Sparsity structure of matrix A']);
    
    figure(2);
    clf;
    spy(Amod,5);
    set(gca,'fontsize',16);
    xlabel('');
    g1=text(0.05*ntmp,0.8*ntmp,['Non-zeros = ',num2str(Amod_nz)]);
    set(g1,'fontsize',16);
    set(g1,'FontWeight','bold');
    title(['Sparsity structure of matrix PAP^{-1}']);
    
    figure(3);
    clf;
    spy(Rmod,5);
    set(gca,'fontsize',16);
    xlabel('');
    g2=text(0.05*ntmp,0.8*ntmp,['Non-zeros = ',num2str(Rmod_nz)]);
    set(g2,'fontsize',16);
    set(g2,'FontWeight','bold');
    title(['Sparsity structure of matrix R']);
    
    figure(4);
    clf;
    spy(Lmod,5);
    set(gca,'fontsize',16);
    xlabel('');
    g2=text(0.05*ntmp,0.8*ntmp,['Non-zeros = ',num2str(Lmod_nz)]);
    set(g2,'fontsize',16);
    set(g2,'FontWeight','bold');
    title(['Sparsity structure of matrix L']);
    
    figure(1);
    
end
