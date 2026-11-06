function probe = fsdt_probe(mesh,point)

point=point(:)';
if numel(point)~=2 || any(~isfinite(point)) || point(1)<0 || ...
        point(1)>mesh.a || point(2)<0 || point(2)>mesh.b
    error('Output point must lie on the plate.');
end
dx=mesh.a/mesh.nx; dy=mesh.b/mesh.ny;
i=min(mesh.nx,max(1,floor(point(1)/dx)+1));
j=min(mesh.ny,max(1,floor(point(2)/dy)+1));
e=(j-1)*mesh.nx+i;
xi=2*(point(1)-(i-1)*dx)/dx-1;
eta=2*(point(2)-(j-1)*dy)/dy-1;
[N,dnat]=shapeFunctionQ4(xi,eta);
nodes=mesh.elements(e,:); [~,~,dN]=Jacobian(mesh.nodes(nodes,:),dnat);
probe.row=sparse(1,nodes,N',1,mesh.ndof);
probe.nodes=nodes; probe.element=e; probe.N=N; probe.dN=dN;
probe.dof=[nodes nodes+mesh.nn nodes+2*mesh.nn nodes+3*mesh.nn nodes+4*mesh.nn];
probe.point=point; probe.natural=[xi eta];
end
