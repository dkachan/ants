
function [node,element,stressPoints] = meshRectangularRegion(pt1, pt2, pt3, pt4, nnx,nny,elemType)

% obtain number of elements from number of nodes


switch elemType
   
    case 'Q4'           % here we generate the mesh of Q4 elements
        numx = nnx-1;
        numy = nny-1;
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);
        stressPoints=[-1 -1;1 -1;1 1;-1 1];  

    case 'Q9'           % here we generate a mesh of Q9 elements
        numx = (nnx-1)/2;
        numy = (nny-1)/2;
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=2;
        inc_v=2*nnx;
        node_pattern=[ 1 3 2*nnx+3 2*nnx+1 2 nnx+3 2*nnx+2 nnx+1 nnx+2 ];
        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);
        stressPoints=[-1 -1;1 -1;1 1;-1 1; 0 -1; 1 0; 0 1; -1 0; 0 0];  
        
    case 'T3' % and last but not least T3 elements
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        node_pattern1=[ 1 2 nnx+1 ];
        node_pattern2=[ 2 nnx+2 nnx+1 ];
        inc_u=1;
        inc_v=nnx;
        element=[make_elem(node_pattern,numx,numy,inc_u,inc_v);
            make_elem(node_pattern,numx,numy,inc_u,inc_v)];
    otherwise
        error('For now Q4, Q9 and T3 are supported by the mesh generator');
end


end
