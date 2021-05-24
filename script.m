%% Clear workspace
commandwindow; clc; clear all; close all

% Check for GMSH
gmsh_path = '/Volumes/gmsh-4.7.1-MacOSX/Gmsh.app/Contents/MacOS/gmsh';
if ~isfile(gmsh_path)
    warning('Download GMSH 4.7.1 from https://gmsh.info/bin/ and provide path to executable "gmsh".');
end

% save_figs = true;
save_figs = false;

% Test cases
%    t: compute solution at this value of time t
%    w(y): interface perturbation function (interface at x = ell + epsilon*w(y))
%    wd(y): derivative of w(y)
%    epsilon: interface perturbation parameter (interface at x = ell + epsilon*w(y))

% Produces Figure 2(a)(b)(c)
% t = 0.2; fig_label = {'(a)','(b)','(c)'}; cnt = 1; w = @(y) zeros(size(y)); wd = @(y) zeros(size(y)); epsilon = 0.0;

% Produces Figure 2(d)(e)(f)
% t = 0.2; fig_label = {'(d)','(e)','(f)'}; cnt = 2; w = @(y) y; wd = @(y) 1; epsilon = 0.05;

% Produces Figure 2(g)(h)(i)
t = 0.2; fig_label = {'(g)','(h)','(i)'}; cnt = 3; w = @(y) sin(pi*y); wd = @(y) pi*cos(pi*y); epsilon = 0.05;

% Produces Figure 2(j)(k)(l)
% t = 0.2; fig_label = {'(j)','(k)','(l)'}; cnt = 4; w = @(y) sin(7*pi*y); wd = @(y) 7*pi*cos(7*pi*y); epsilon = 0.02;

%% Parameters
L = 1; % domain length
H = 1; % domain width
ell = L/2; % unperturbed interface location
D1 = 1; % diffusivity in region 1 (0 < x < l + epsilon*w(y))
D2 = 0.01; % diffusivity in region 2 (l + epsilon*w(y) < x < L)
C0 = @(s) 1./s; % boundary condition function at x = 0
CL = @(s) 0.*s; % boundary condition function at x = L
N = 5; % number of terms in perturbation solution
M = 30; % number of terms in each eigenfunction expansions
ref = 0.01; % refinement parameter in meshing (controls number of elements/nodes)
intrfcN = 101; % number of nodes to use on the finite volume mesh interface
Np = 14; % number of poles/residues used in the numerical inverse Laplace transform

% Check for presence of colormaps
if ~isfile('plasma.m') || ~isfile('viridis.m')
    warning(['Download plasma and viridis colormaps from https://www.mathworks.com/matlabcentral/fileexchange',...
        '/62729-matplotlib-perceptually-uniform-colormaps and place in current directory.'])
    solncmap = parula(128); % colormap for solutions
    errcmap = summer(128); % colormap for difference    
else
    solncmap = plasma(128); % colormap for solutions
    errcmap = viridis(128); % colormap for difference
end

%% Finite Volume Method Solution
% Meshing
generate_mesh(L,H,ell,epsilon,w,ref,intrfcN,gmsh_path);
mesh;
nodes = msh.POS(:,1:2);
elements = msh.TRIANGLES(:,1:3);
element_material = msh.TRIANGLES(:,4);
dirichlet_nodes = msh.LINES(:,1:2);
dirichlet_nodes = unique(dirichlet_nodes(:));
mesh = mesh_properties(nodes,elements,dirichlet_nodes);
D = element_material*D1 + (1-element_material)*D2; % diffusivity per element

% Compute solution
u0 = zeros(size(nodes,1),1); % initial condition
un = fvmsol(D,mesh,t,u0); 

%% Perturbation Solution
[z,c] = cf(Np); % inverse Laplace transform weights/poles

idx = nodes(:,1) <= ell + epsilon*w(nodes(:,2)); % interface nodes
nodes1 = nodes(idx,:); % nodes in region 1 (including interface)
nodes2 = nodes(~idx,:); % nodes in region 2

% Leading order term
[U1,U2] = leading_order_term(D1,D2,ell,L,C0,CL);
u1 = @(t) inverse_laplace_transform(@(s) U1(nodes1(:,1),s,0),t,Np,c,z);
u2 = @(t) inverse_laplace_transform(@(s) U2(nodes2(:,1),s,0),t,Np,c,z);
u0 = @(t) [u1(t); u2(t)];

% Higher order terms
tfun = @(s) higher_order_terms(N,w,wd,U1,U2,D1,D2,H,L,ell,nodes1(:,1),nodes2(:,1),...
    nodes1(:,2),nodes2(:,2),M,s);
ui = inverse_laplace_transform(tfun,t,Np,c,z);

% Asymptotic expansions [Equations (5) and (6)]
u = zeros(size(nodes,1),1); % Pre-allocate the perturbation solution
upert = u0(t);
for i = 1:N-1
    upert = upert + epsilon^i*ui(:,i);
end
u(idx) = upert(1:sum(idx));
u(~idx) = upert(sum(idx)+1:end);

%% Comparison
figure
set(gcf,'position',[200,200,560*3.1 420],'Color','w')

% Perturbation solution
p1 = subplot(1,3,1);
trisurf(elements,nodes(:,2),L-nodes(:,1),u);
colormap(p1,solncmap)
shading interp
hold on
zlim([0,1])
caxis([0,1])
set(gca,'DataAspectRatio',[1,1,1],'FontSize',24,'Xtick',[],'Ytick',[])
text(-0.07,-0.07,fig_label{1},'FontSize',24)
view(2)
colorbar('Ticks',[0:0.1:1]);
box on

% Finite Volume Method solution
p2 = subplot(1,3,2);
trisurf(elements,nodes(:,2),L-nodes(:,1),un,'LineStyle','-');
colormap(p2,solncmap)
shading interp
hold on
set(gca,'DataAspectRatio',[1,1,1],'FontSize',24,'Xtick',[],'Ytick',[])
text(-0.07,-0.07,fig_label{2},'FontSize',24)
view(2)
colorbar('Ticks',[0:0.1:1]);
text(1.115,0,'0','FontSize',21)
box on

% Difference between the solutions
p3 = subplot(1,3,3);
trisurf(elements,nodes(:,2),L-nodes(:,1),log10(abs(u-un)));
colormap(p3,errcmap)
shading interp
zlim([-16,-2])
caxis([-16,-2])
set(gca,'DataAspectRatio',[1,1,1],'FontSize',24,'Xtick',[],'Ytick',[])
text(-0.07,-0.07,fig_label{3},'FontSize',24)
view(2)
colorbar('Ticks',[-16:2:-2],'YTickLabel',{'10^{-16}','10^{-14}','10^{-12}','10^{-10}',...
    '10^{-8}','10^{-6}','10^{-4}','10^{-2}'});
box on
fprintf('\n%1.2e\n',max(max(abs(u-un))))

%% Save figs
% addpath(insertpath)
% path_name = '../../Figures/';
% name = 'Case';
% if save_figs
%     feval('export_fig',gcf,[path_name,name,num2str(cnt)],'-pdf')
% end
