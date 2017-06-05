function A = ConstructA_CG3(SubNumPhysNodes,NumPhysElems,...
                            SubNumBndNodes,sub_node,...
                            tnode,node_subs,area,sub_bnd_node,Jmat)
  
    A = sparse(SubNumPhysNodes,SubNumPhysNodes);
    
    A1 = zeros(10,1);
    A2 = zeros(10,1);
    A3 = zeros(10,1);
    A4 = zeros(10,1);
    A5 = zeros(10,1);
    A6 = zeros(10,1);
    A7 = zeros(10,1);
    A8 = zeros(10,1);
    A9 = zeros(10,1);
    A10 = zeros(10,1);
    
    A1(1) =   0.0;
    A1(2) =   0.0;
    A1(3) =   0.0;
    A1(4) =   0.0;
    A1(5) =   0.0;
    A1(6) =   1.0;
    A1(7) =   0.0;
    A1(8) =   0.0;
    A1(9) =   0.0;
    A1(10) =  0.0;
    
    A2(1) =   0.5;
    A2(2) =  -1.5;
    A2(3) =   1.5;
    A2(4) =  -0.5;
    A2(5) =  -1.5;
    A2(6) =   0.0;
    A2(7) =   1.5;
    A2(8) =   0.0;
    A2(9) =   0.0;
    A2(10) =  0.0;
    
    A3(1) =   0.5;
    A3(2) =  -1.5;
    A3(3) =   0.0;
    A3(4) =   0.0;
    A3(5) =  -1.5;
    A3(6) =   0.0;
    A3(7) =   0.0;
    A3(8) =   1.5;
    A3(9) =   1.5;
    A3(10) = -0.5;
    
    A4(1) =   0.0;
    A4(2) =   4.5;
    A4(3) =  -4.5;
    A4(4) =   0.0;
    A4(5) =   4.5;
    A4(6) =  -9.0;
    A4(7) =   4.5;
    A4(8) =  -4.5;
    A4(9) =   4.5;
    A4(10) =  0.0;
    
    A5(1) =   0.0;
    A5(2) =   0.0;
    A5(3) =   0.0;
    A5(4) =   0.0;
    A5(5) =   4.5;
    A5(6) =  -9.0;
    A5(7) =   4.5;
    A5(8) =   0.0;
    A5(9) =   0.0;
    A5(10) =  0.0;
    
    A6(1) =   0.0;
    A6(2) =   4.5;
    A6(3) =   0.0;
    A6(4) =   0.0;
    A6(5) =   0.0;
    A6(6) =  -9.0;
    A6(7) =   0.0;
    A6(8) =   0.0;
    A6(9) =   4.5;
    A6(10) =  0.0;
  
    A7(1) =  -4.5;
    A7(2) =  13.5;
    A7(3) = -13.5;
    A7(4) =   4.5;
    A7(5) =   0.0;
    A7(6) =   0.0;
    A7(7) =   0.0;
    A7(8) =   0.0;
    A7(9) =   0.0;
    A7(10) =  0.0;
    
    A8(1) = -13.5;
    A8(2) =  27.0;
    A8(3) = -13.5;
    A8(4) =   0.0;
    A8(5) =  13.5;
    A8(6) = -27.0;
    A8(7) =  13.5;
    A8(8) =   0.0;
    A8(9) =   0.0;
    A8(10) =  0.0;
    
    A9(1) = -13.5;
    A9(2) =  13.5;
    A9(3) =   0.0;
    A9(4) =   0.0;
    A9(5) =  27.0;
    A9(6) = -27.0;
    A9(7) =   0.0;
    A9(8) = -13.5;
    A9(9) =  13.5;
    A9(10) =  0.0;
    
    A10(1) =  -4.5;
    A10(2) =   0.0;
    A10(3) =   0.0;
    A10(4) =   0.0;
    A10(5) =  13.5;
    A10(6) =   0.0;
    A10(7) =   0.0;
    A10(8) = -13.5;
    A10(9) =   0.0;
    A10(10) =  4.5;
    
    spts = zeros(6,2);
    spts(1,1) =  0.112615157582632;
    spts(1,2) =  0.112615157582632;
    
    spts(2,1) = -0.225230315165263;
    spts(2,2) =  0.112615157582632;
    
    spts(3,1) =  0.112615157582632;
    spts(3,2) = -0.225230315165263;
    
    spts(4,1) = -0.241757119823562;
    spts(4,2) = -0.241757119823562;
    
    spts(5,1) =  0.483514239647126;
    spts(5,2) = -0.241757119823562;
    
    spts(6,1) = -0.241757119823562;
    spts(6,2) =  0.483514239647126;
    
    wgts = zeros(6,1);
    wgts(1) = 0.1116907948390055;
    wgts(2) = 0.1116907948390055;
    wgts(3) = 0.1116907948390055;
    wgts(4) = 0.0549758718276610;
    wgts(5) = 0.0549758718276610;
    wgts(6) = 0.0549758718276610;
    
    % Loop over all elements in the mesh
    for i=1:NumPhysElems
    
        % Information for element i
        tt = node_subs(i,:);
                    
        % Evaluate gradients of the Lagrange polynomials on Gauss quadrature points      
        gpx = zeros(10,6);
        gpy = zeros(10,6);
        
        for m=1:6

            xi  = spts(m,1);
            eta = spts(m,2);
            
            for k=1:10
                
                gp_xi  = A2(k) + A4(k)*eta + 2.0*A5(k)*xi ...
                         + 3.0*A7(k)*xi*xi + 2.0*A8(k)*xi*eta + A9(k)*eta*eta;
                
                gp_eta = A3(k) + A4(k)*xi + 2.0*A6(k)*eta ...
                         + A8(k)*xi*xi + 2.0*A9(k)*xi*eta + 3.0*A10(k)*eta*eta;

                gpx(k,m) = Jmat(i,1,1)*gp_xi + Jmat(i,1,2)*gp_eta;
                gpy(k,m) = Jmat(i,2,1)*gp_xi + Jmat(i,2,2)*gp_eta;

            end
        end
        
        % Entries of the stiffness matrix A
        T = area(i);
        for j=1:10
            for k=1:10
                
                tmp = A(tt(j),tt(k));
                
                for m=1:6
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