function A = ConstructA_CG1(NumPhysNodes,NumPhysElems,NumBndNodes,node,tnode,area,bnd_node)

    A = sparse(NumPhysNodes,NumPhysNodes);

    % Loop over all elements in the mesh
    for i=1:NumPhysElems

        % Information for element i
        tt = tnode(i,:);
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

                tmp = A(tt(j),tt(k)) + T*(gphi(j,1)*gphi(k,1)+ ...
                                          gphi(j,2)*gphi(k,2));
                A(tt(j),tt(k)) = tmp;

            end
        end
    end

    % Replace boundary node equations by Dirichlet boundary condition enforcement
    for i=1:NumBndNodes
        
        j = bnd_node(i);

        A(j,:) = 0;
        A(:,j) = 0;
        A(j,j) = 1;

    end