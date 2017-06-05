%
%  Add extra points and elements if points_per_dir>1
%
function [node,tnode,z]=DivideUnst2Mesh(points_per_dir,NumPhysElems,node,tnode)
%DIVIDEUNST2MESH.    Divide the unstructured grid into a finer mesh.
%
% This function is used to plot 'multiple' points per element.
%
% Create a new mesh with higher resolution.  "pretends that it is piecewise
% constant on each 'element', which is actually a subdivision of each element.
  
    k  = 0;  % new node index
    kk = 0;
    m  = 0;  % new element index
    
    node_list = [];
    map_list  = [];

    for n=1:NumPhysElems
      
      h = 1/points_per_dir;
      
      x1 = node(tnode(n,1),1);
      y1 = node(tnode(n,1),2);
      
      x2 = node(tnode(n,2),1);
      y2 = node(tnode(n,2),2);
      
      x3 = node(tnode(n,3),1);
      y3 = node(tnode(n,3),2);
      
      xm = (x1+x2+x3)/3;
      ym = (y1+y2+y3)/3;
      
      km = zeros((points_per_dir+1)*(points_per_dir+2)/2,1);
      index = 0;
      
      for j=1:(points_per_dir+1)
        for i=1:(points_per_dir+2-j) 
          index = index+1;      
          k = k+1;
          km(index) = k;
          
          xi  = (i-1)*h-1/3;
          eta = (j-1)*h-1/3;
          
          xtmp = xm + (x2-x1)*xi + (x3-x1)*eta;
          ytmp = ym + (y2-y1)*xi + (y3-y1)*eta;
          
          ztmp = 3*pi + atan2(ytmp,xtmp) + 1.0e3*sqrt(xtmp^2+ytmp^2);
          
          if (isempty(node_list)==1)
            kk = 1;
            node_list = ztmp;
            map_list = kk;
            
            node_new(kk,1) = xtmp;
            node_new(kk,2) = ytmp;     
          else
            [ktmp,ltmp] = min(abs(ztmp-node_list));
            if ( ktmp >= 1.0e-8 )
              kk = kk+1;
              node_list = [node_list,ztmp];
              map_list  = [map_list,kk];
              node_new(kk,1) = xtmp;
              node_new(kk,2) = ytmp;
            else
              map_list  = [map_list,ltmp];                        
            end
          end
        end
      end
      
      % all triangles oriented one way     
      m1 = km(1);
      m2 = km(points_per_dir+2);
      for j=1:points_per_dir
        for i=1:(points_per_dir+1-j)
          m = m+1;            
          tnode_new(m,1) = map_list(m1+i-1);
          tnode_new(m,2) = map_list(m1+i);
          tnode_new(m,3) = map_list(m2+i-1);
        end
        m1 = m1 + (points_per_dir+2-j);
        m2 = m2 + (points_per_dir+1-j);
      end 
      % all triangles oriented the other way
      m1 = km(1);
      m2 = km(points_per_dir+2);
      for j=1:(points_per_dir-1)
        for i=1:(points_per_dir-j)
          m = m+1;            
          tnode_new(m,1) = map_list(m1+i);
          tnode_new(m,2) = map_list(m2+i);
          tnode_new(m,3) = map_list(m2+i-1);
        end       
        m1 = m1 + (points_per_dir+2-j);
        m2 = m2 + (points_per_dir+1-j);
      end
    
    end
                
    % Create element midpoints in canonical coordinates
    z = zeros(points_per_dir^2, 2);     
    
    x1 = node(tnode(1,1),1);
    y1 = node(tnode(1,1),2);
    
    x2 = node(tnode(1,2),1);
    y2 = node(tnode(1,2),2);
    
    x3 = node(tnode(1,3),1);
    y3 = node(tnode(1,3),2);

    Area = (x3-x2)*y1 + (x1-x3)*y2 + (x2-x1)*y3;    
    
    for m=1:(points_per_dir^2)
      
      p1 = (node_new(tnode_new(m,1),1)+node_new(tnode_new(m,2),1)+...
            node_new(tnode_new(m,3),1))/3;
      
      p2 = (node_new(tnode_new(m,1),2)+node_new(tnode_new(m,2),2)+...
            node_new(tnode_new(m,3),2))/3;
      
      z(m,1) = ((x1-x3)*p2+(-y1+y3)*p1+(2/3)*y1*x3+(1/3)*y1*x2-...
               (1/3)*y2*x1-(1/3)*y3*x2-(2/3)*y3*x1+(1/3)*y2*x3)/Area;
      
      z(m,2) = ((-x1+x2)*p2+(y1-y2)*p1+(2/3)*y2*x1+(1/3)*y3*x1-...
               (2/3)*y1*x2-(1/3)*y1*x3+(1/3)*y2*x3-(1/3)*y3*x2)/Area;
    
    end      
    
    tnode = tnode_new;
    node = node_new;
    clear tnode_new;
    clear node_new;
           
    % Properly orient each triangle    
    for j=1:size(tnode,1)
    
      x1 = node(tnode(j,1),1);
      y1 = node(tnode(j,1),2);
      
      x2 = node(tnode(j,2),1);
      y2 = node(tnode(j,2),2);
      
      x3 = node(tnode(j,3),1);
      y3 = node(tnode(j,3),2);

      Area = (x3-x2)*y1 + (x1-x3)*y2 + (x2-x1)*y3;
    
      if (Area<0)
        tmp = tnode(j,2);
        tnode(j,2) = tnode(j,3);
        tnode(j,3) = tmp;
      end
      
    end

end
