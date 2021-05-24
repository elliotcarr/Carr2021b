function generate_mesh(L,H,ell,epsilon,w,ref,N,gmsh_path)

% Define geometry in gmsh
fid = fopen('mesh.geo', 'w');
fprintf(fid,'r = %g;\n',ref);
fprintf(fid,'Mesh.Algorithm = 6;\n');
fprintf(fid,'Mesh.Format = 50;\n'); % matlab output (creates mesh.m)
fprintf(fid,'Point(1) = {%g,%g,0,r};\n',0,0);
fprintf(fid,'Point(2) = {%g,%g,0,r};\n',L,0);
fprintf(fid,'Point(3) = {%g,%g,0,r};\n',L,H);
fprintf(fid,'Point(4) = {%g,%g,0,r};\n',0,H);
y = linspace(0,H,N);
x = ell + epsilon*w(y);
for ii = 1:N
    fprintf(fid,'Point(%i) = {%g,%g,0,r};\n',ii+4,x(ii),y(ii));
end
fprintf(fid,'Line(1) = {1,5};\n');
fprintf(fid,'Line(2) = {5,2};\n');
fprintf(fid,'Line(3) = {2,3};\n');
fprintf(fid,'Line(4) = {3,%i};\n',N+4);
fprintf(fid,'Line(5) = {%i,4};\n',N+4);
fprintf(fid,'Line(6) = {4,1};\n');
fprintf(fid,'BSpline(7) = {');
for ii = (1:N)+4
    if ii == N+4
        fprintf(fid,'%i',ii);
    else
        fprintf(fid,'%i,',ii);
    end
end
fprintf(fid,'};\n');
fprintf(fid,'Line Loop(1) = {1,7,5,6};\n');
fprintf(fid,'Line Loop(2) = {2,3,4,-7};\n');
fprintf(fid,'Physical Line(99) = {6};\n'); % label nodes along line 6 (Dirichlet nodes)
fprintf(fid,'Plane Surface(1) = {1};\n');
fprintf(fid,'Plane Surface(2) = {2};\n');
fprintf(fid,'Physical Surface(1) = {1};\n');
fprintf(fid,'Physical Surface(0) = {2};\n');
fclose(fid);

% Performing meshing
system([gmsh_path,' mesh.geo -2']);

end
