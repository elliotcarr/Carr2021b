function mesh = mesh_properties(nodes,elements,dirichlet_nodes)

no_elements  = size(elements,1);
no_nodes     = size(nodes,1);

% Initialise
V          = zeros(no_nodes,1);
nx         = zeros(no_elements,3);
ny         = zeros(no_elements,3);
s          = zeros(no_elements,6);
L          = zeros(no_elements,3);

% Compute geometrical properties of mesh
for k = 1:no_elements
    
    v = elements(k,:); % vertices of kth element
    xv   = nodes(v,1); % x-coordinates of the three nodes/vertices
    yv   = nodes(v,2); % y-coordinates of the three nodes/vertices
    
    % Element centroid
    centroid = [sum(xv)/3,sum(yv)/3];
    
    % Element edge midpoint
    m1   = [(xv(1)+xv(2))/2,(yv(1)+yv(2))/2];
    m2   = [(xv(2)+xv(3))/2,(yv(2)+yv(3))/2];
    m3   = [(xv(3)+xv(1))/2,(yv(3)+yv(1))/2];
    
    % Control volume edge vectors
    e1 = centroid - m1;
    e2 = centroid - m2;
    e3 = centroid - m3;
    
    % Control volume edge lengths
    L(k,1) = norm(e1,2);
    L(k,2) = norm(e2,2);
    L(k,3) = norm(e3,2);

    % Unit normals
    nx(k,1) = e1(2)/L(k,1); ny(k,1) = -e1(1)/L(k,1);
    nx(k,2) = e2(2)/L(k,2); ny(k,2) = -e2(1)/L(k,2);
    nx(k,3) = e3(2)/L(k,3); ny(k,3) = -e3(1)/L(k,3);
    
    % Vectors connecting nodes to centroid
    pj = centroid - [xv,yv];
    p1 = pj(1,:);
    p2 = pj(2,:);
    p3 = pj(3,:);
    
    % Vectors connecting midpoints
    q1 = m1 - m3;
    q2 = m2 - m1;
    q3 = m3 - m2;
    
    % Sub-control volume areas
    S1 = 0.5*(q1(1)*p1(2) - q1(2)*p1(1));
    S2 = 0.5*(q2(1)*p2(2) - q2(2)*p2(1));
    S3 = 0.5*(q3(1)*p3(2) - q3(2)*p3(1));
    
    % Control volume areas
    V(v(1)) = V(v(1)) + S1;
    V(v(2)) = V(v(2)) + S2;
    V(v(3)) = V(v(3)) + S3;

    % Shape function information (linear interpolation)
    detA = xv(2)*yv(3) - yv(2)*xv(3) - xv(1)*yv(3) + xv(3)*yv(1) + xv(1)*yv(2) - xv(2)*yv(1);
    s(k,1) = (yv(2)-yv(3)) / detA;
    s(k,2) = (yv(3)-yv(1)) / detA;
    s(k,3) = (yv(1)-yv(2)) / detA;
    s(k,4) = (xv(3)-xv(2)) / detA;
    s(k,5) = (xv(1)-xv(3)) / detA;
    s(k,6) = (xv(2)-xv(1)) / detA;
    
end

mesh.elements        = elements;
mesh.nodes           = nodes;
mesh.no_elements     = no_elements;
mesh.no_nodes        = no_nodes;
mesh.nx              = nx;
mesh.ny              = ny;
mesh.s               = s;
mesh.V               = V;
mesh.L               = L;
mesh.dirichlet_nodes = dirichlet_nodes;

end
