function mesh = fsdt_mesh(a,b,nx,ny)
%FSDT_MESH Open replacement for Ferreira's supplied rectangularMesh.p.
% Counterclockwise Q4: bottom-left, bottom-right, top-right, top-left.
if any(~isfinite([a b nx ny])) || a<=0 || b<=0 || ...
        nx<2 || ny<2 || nx~=round(nx) || ny~=round(ny)
    error('Positive plate dimensions and integer nx,ny >= 2 required.');
end
[X,Y]=ndgrid(linspace(0,a,nx+1),linspace(0,b,ny+1));
mesh.nodes=[X(:) Y(:)]; mesh.elements=zeros(nx*ny,4);
e=0;
for j=1:ny
    for i=1:nx
        e=e+1; n=(j-1)*(nx+1)+i;
        mesh.elements(e,:)=[n n+1 n+nx+2 n+nx+1];
    end
end
mesh.a=a; mesh.b=b; mesh.nx=nx; mesh.ny=ny;
mesh.nn=size(mesh.nodes,1); mesh.ne=size(mesh.elements,1);
mesh.ndof=5*mesh.nn;
end
