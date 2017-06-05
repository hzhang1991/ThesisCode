function A = ConstructAperiodic_CG1(NumPhysNodes,...
                                    NumPhysElems,...
                                    NumBndNodes,...
                                    NumRdntBndNodes,...
                                    node,...
                                    tnode,...
                                    area,...
                                    bnd_node,...
                                    node_map,...
                                    rdnt_bnd_node)

    A = sparse(NumPhysNodes,NumPhysNodes);

    % Loop over all elements in the mesh
    for i=1:NumPhysElems

        % Information for element i
        tt = tnode(i,:);
        ss = node_map(tt);
        x = node(tt,1);
        y = node(tt,2);
        T = area(i);

        % Compute the three Lagrange polynomials on this element
        Mat = [1,x(1),y(1);
               1,x(2),y(2);
               1,x(3),y(3)];
        MatInv = inv(Mat);

        % Compute the gradients of the three Lagrange polynomials
        % on this element
        gphi = sparse(3,2);
        gphi(1:3,1) = transpose(MatInv(2,1:3));
        gphi(1:3,2) = transpose(MatInv(3,1:3));

        % Entries of the stiffness matrix A
        for j=1:3
            for k=1:3

                tmp = A(ss(j),ss(k)) + T*(gphi(j,1)*gphi(k,1)+ ...
                                          gphi(j,2)*gphi(k,2));
                A(ss(j),ss(k)) = tmp;

            end
        end
    end
    
    for i=1:NumRdntBndNodes
        
        j = rdnt_bnd_node(i,1);
        A(j,j) = 1;
        
    end
    
    A(NumPhysNodes,NumPhysNodes) = 1+A(NumPhysNodes,NumPhysNodes);