function un = fvmsol(D,mesh,t,u0)
% Computes finite volume solution at time t
% Note this code only works for the boundary conditions c0(t) = 1 and qL(t) = 0.

% Define relevant properties from the mesh
elements        = mesh.elements;
no_elements     = mesh.no_elements;
no_nodes        = mesh.no_nodes;
nx              = mesh.nx;
ny              = mesh.ny;
s               = mesh.s;
V               = mesh.V;
L               = mesh.L;
dirichlet_nodes = mesh.dirichlet_nodes;

% Build A and b for ODE system: du/dt = Au + b
A = zeros(no_nodes);
b = zeros(no_nodes,1);
for k = 1:no_elements
    
    % Vertices of kth element
    v = elements(k,:);
    
    for j = 1:3
        
        if j == 1
            jnb = 2;
        elseif j == 2
            jnb = 3;
        elseif j == 3
            jnb = 1;
        end
        
        % Linear interpolation for gradient approximation
        Av1 = -D(k)*(s(k,1)*nx(k,j) + s(k,4)*ny(k,j))*L(k,j);
        Av2 = -D(k)*(s(k,2)*nx(k,j) + s(k,5)*ny(k,j))*L(k,j);
        Av3 = -D(k)*(s(k,3)*nx(k,j) + s(k,6)*ny(k,j))*L(k,j);
        
        A(v(j),v(1)) = A(v(j),v(1)) - Av1;
        A(v(j),v(2)) = A(v(j),v(2)) - Av2;
        A(v(j),v(3)) = A(v(j),v(3)) - Av3;
        
        A(v(jnb),v(1)) = A(v(jnb),v(1)) + Av1;
        A(v(jnb),v(2)) = A(v(jnb),v(2)) + Av2;
        A(v(jnb),v(3)) = A(v(jnb),v(3)) + Av3;
        
    end
    
end

% Scale each node by the size of its respective control volume area
A = A./V;

% Remove nodes on Dirichlet boundary (x=0) from ODE system 
% (as solution is known)
for i = 1:length(dirichlet_nodes)
    b = b + A(:,dirichlet_nodes(i))*1;
end
A(dirichlet_nodes,:) = [];
A(:,dirichlet_nodes) = [];
b(dirichlet_nodes) = [];
u0(dirichlet_nodes) = [];

% Solve initial value problem using ode15s
A = sparse(A);
b = sparse(b);
options = odeset('Jacobian',A,'Jpattern',A,'JConstant',true,'RelTol',1e-12,'AbsTol',1e-12);
[~,ut] = ode15s(@(t,u) A*u + b,[0,t],u0,options);
ut = ut';
ut(:,1:end-1) = [];

% Add back in solution on Dirichlet boundary (x = 0)
var_num = zeros(no_nodes,1);
cnt = 0;
for i = 1:no_nodes
    if ismember(i,dirichlet_nodes)
        var_num(i) = 0;
    else 
        cnt = cnt + 1;
        var_num(i) = cnt;  
    end
end
un = zeros(no_nodes,1);
for i = 1:no_nodes
    if var_num(i) == 0
        un(i,:) = 1;
    else
        un(i,:) = ut(var_num(i),:);
    end
end
