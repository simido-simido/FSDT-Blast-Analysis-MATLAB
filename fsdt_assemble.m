function model = fsdt_assemble(cfg)

mesh=fsdt_mesh(cfg.plate.a,cfg.plate.b,cfg.mesh.nx,cfg.mesh.ny);
lam=fsdt_laminate(cfg.plies,cfg.shearCorrection);
patch=cfg.load.patch;
if isempty(patch), patch=[0 mesh.a 0 mesh.b]; end
if numel(patch)~=4 || any(~isfinite(patch)) || ...
        patch(1)<0 || patch(2)>mesh.a || patch(3)<0 || patch(4)>mesh.b || ...
        patch(2)<=patch(1) || patch(4)<=patch(3)
    error('load.patch must be [] or [xmin xmax ymin ymax] inside plate.');
end
tol=1e-9;
ratios=[patch(1:2)/(mesh.a/mesh.nx) patch(3:4)/(mesh.b/mesh.ny)];
if any(abs(ratios-round(ratios))>tol)
    error('Pressure patch must follow element boundaries. Refine/align mesh.');
end
if ~isscalar(cfg.load.direction) || ~ismember(cfg.load.direction,[-1 1])
    error('load.direction must be +1 or -1.');
end
entries=400*mesh.ne; ii=zeros(entries,1); jj=ii; kk=ii; mm=ii;
fp=zeros(mesh.ndof,1);
[weights,points]=gaussQuadrature('complete');
for e=1:mesh.ne
    nodes=mesh.elements(e,:); xy=mesh.nodes(nodes,:);
    dof=[nodes nodes+mesh.nn nodes+2*mesh.nn nodes+3*mesh.nn nodes+4*mesh.nn];
    Ke=zeros(20); Me=zeros(20); fe=zeros(4,1);
    center=mean(xy,1);
    loaded=center(1)>patch(1) && center(1)<patch(2) && ...
        center(2)>patch(3) && center(2)<patch(4);
    for g=1:4
        [N,dnat]=shapeFunctionQ4(points(g,1),points(g,2));
        [J,~,dN]=Jacobian(xy,dnat);
        detJ=det(J); if detJ<=0, error('Non-positive element Jacobian.'); end
        [Bm,Bb,~]=fsdt_B(N,dN); dv=weights(g)*detJ;
        Ke=Ke+(Bm'*lam.A*Bm+Bm'*lam.B*Bb+Bb'*lam.B'*Bm+Bb'*lam.D*Bb)*dv;
        Me=Me+kron(lam.inertia,N*N')*dv;
        if loaded, fe=fe+N*dv; end
    end
    [N,dnat]=shapeFunctionQ4(0,0); [J,~,dN]=Jacobian(xy,dnat);
    [~,~,Bs]=fsdt_B(N,dN); Ke=Ke+Bs'*lam.S*Bs*4*det(J);
    [rr,cc]=ndgrid(dof,dof); ix=(e-1)*400+(1:400);
    ii(ix)=rr(:); jj(ix)=cc(:); kk(ix)=Ke(:); mm(ix)=Me(:);
    fp(nodes)=fp(nodes)+fe;
end
model.K=sparse(ii,jj,kk,mesh.ndof,mesh.ndof);
model.M=sparse(ii,jj,mm,mesh.ndof,mesh.ndof);
model.forceUnit=cfg.load.direction*fp;
model.mesh=mesh; model.lam=lam; model.patch=patch;
[model.fixed,model.free]=fsdt_boundary(mesh,cfg.boundary);
model.Kff=model.K(model.free,model.free); model.Mff=model.M(model.free,model.free);
model.forceUnitFree=model.forceUnit(model.free);
model.loadedArea=(patch(2)-patch(1))*(patch(4)-patch(3));
if abs(sum(fp)-model.loadedArea)>1e-10*model.loadedArea
    error('Pressure force does not integrate to loaded area.');
end
end
