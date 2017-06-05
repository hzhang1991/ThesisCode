function A = ConstructA_CG2(SubNumPhysNodes,NumPhysElems,...
                            SubNumBndNodes,sub_node,...
                            tnode,node_subs,area,sub_bnd_node,Jmat)
  
    A = sparse(SubNumPhysNodes,SubNumPhysNodes);
    
    A1 = zeros(6,1);
    A2 = zeros(6,1);
    A3 = zeros(6,1);
    A4 = zeros(6,1);
    A5 = zeros(6,1);
    A6 = zeros(6,1);
    
    oneninth = 1/9;
    onethird = 1/3;
    
    A1(1) = -oneninth;
    A1(2) =  4.0*oneninth;
    A1(3) = -oneninth;
    A1(4) =  4.0*oneninth;
    A1(5) =  4.0*oneninth;
    A1(6) = -oneninth;
  
    A2(1) = -onethird;
    A2(2) =  0.0;
    A2(3) =  onethird;
    A2(4) = -4.0*onethird;
    A2(5) =  4.0*onethird;
    A2(6) =  0.0;
  
    A3(1) = -onethird;
    A3(2) = -4.0*onethird;
    A3(3) =  0.0;
    A3(4) =  0.0;
    A3(5) =  4.0*onethird;
    A3(6) =  onethird;
    
    A4(1) =  4.0;
    A4(2) = -4.0;
    A4(3) =  0.0;
    A4(4) = -4.0;
    A4(5) =  4.0;
    A4(6) =  0.0;
    
    A5(1) =  2.0;
    A5(2) = -4.0;
    A5(3) =  2.0;
    A5(4) =  0.0;
    A5(5) =  0.0;
    A5(6) =  0.0;
    
    A6(1) =  2.0;
    A6(2) =  0.0;
    A6(3) =  0.0;
    A6(4) = -4.0;
    A6(5) =  0.0;
    A6(6) =  2.0;

    spts = zeros(3,2);
    spts(1,1) =  1.0/3.0;
    spts(1,2) = -1.0/6.0;
    spts(2,1) = -1.0/6.0;
    spts(2,2) = -1.0/6.0;
    spts(3,1) = -1.0/6.0;
    spts(3,2) =  1.0/3.0;
  
    wgts = zeros(3,1);
  
    wgts(1) = 1.0/6.0;
    wgts(2) = 1.0/6.0;
    wgts(3) = 1.0/6.0;
  
    % Loop over all elements in the mesh
    for i=1:NumPhysElems
    
        % Information for element i
        tt = node_subs(i,:);
                    
        % Evaluate gradients of the Lagrange polynomials on Gauss quadrature points      
        gpx = zeros(6,3);
        gpy = zeros(6,3);

        for m=1:3

            xi  = spts(m,1);
            eta = spts(m,2);
            
            for k=1:6

                gp_xi  = A2(k) + 2.0*A5(k)*xi + A4(k)*eta;
                gp_eta = A3(k) + A4(k)*xi + 2.0*A6(k)*eta;

                gpx(k,m) = Jmat(i,1,1)*gp_xi + Jmat(i,1,2)*gp_eta;
                gpy(k,m) = Jmat(i,2,1)*gp_xi + Jmat(i,2,2)*gp_eta;

            end
        end
        
        % Entries of the stiffness matrix A
        T = area(i);
        for j=1:6
            for k=1:6
                
                tmp = A(tt(j),tt(k));
                
                for m=1:3
                    tmp = tmp + 2.0*T*wgts(m)*(gpx(j,m)*gpx(k,m)+gpy(j,m)*gpy(k,m));
                end
                
                A(tt(j),tt(k)) = tmp;
                
            end
        end
        
    end
                

    % Replace boundary node equations by Dirichlet boundary condition enforcement
    for i=1:SubNumBndNodes
        
        j = sub_bnd_node(i);

        A(j,:) = 0;
        A(:,j) = 0;
        A(j,j) = 1;

    end

