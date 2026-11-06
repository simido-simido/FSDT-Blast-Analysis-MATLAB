function modal = fsdt_modes(model,count)
K=model.Kff; M=model.Mff; nf=size(K,1);
if count<1 || count~=round(count), error('Positive integer modal count required.'); end
nreq=min(nf-1,max(20,3*count));
if nf<=350
    [V,D]=eig(full(K),full(M)); vals=real(diag(D));
else
    opts.issym=true; opts.isreal=true; opts.tol=1e-9; opts.maxit=3000;
    opts.v0=sin((1:nf)'*sqrt(2));
    [V,D,flag]=eigs(K,M,nreq,0,opts); vals=real(diag(D));
    if flag~=0, error('Sparse eigensolver did not converge.'); end
end
[vals,ix]=sort(vals); V=real(V(:,ix));
if any(~isfinite(vals)) || vals(1)<=0
    error('Zero/negative eigenvalue: inspect supports and reduced-integration mechanisms.');
end
nkeep=min(count,numel(vals)); vals=vals(1:nkeep); V=V(:,1:nkeep);
for j=1:nkeep, V(:,j)=V(:,j)/sqrt(V(:,j)'*M*V(:,j)); end
allV=zeros(model.mesh.ndof,nkeep); allV(model.free,:)=V;
frac=zeros(nkeep,1);
for j=1:nkeep
    w=allV(1:model.mesh.nn,j);
    frac(j)=real(w'*model.M(1:model.mesh.nn,1:model.mesh.nn)*w);
end
modal.omega=sqrt(vals); modal.frequencyHz=modal.omega/(2*pi);
modal.vectors=allV; modal.vectorsFree=V;
modal.transverseMassFraction=frac;
modal.flexural=find(frac>0.10);
modal.residual=zeros(nkeep,1);
for j=1:nkeep
    kv=K*V(:,j); mv=vals(j)*M*V(:,j);
    modal.residual(j)=norm(kv-mv)/(norm(kv)+norm(mv));
end
end
