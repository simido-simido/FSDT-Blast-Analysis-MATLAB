function ref = fsdt_navier_reference(lam,a,b,m,n)
if norm(lam.B,'fro')>1e-10*sqrt(norm(lam.A,'fro')*norm(lam.D,'fro')) || ...
        norm(lam.D([1 2],3))>1e-10*norm(lam.D,'fro') || ...
        abs(lam.S(1,2))>1e-10*norm(lam.S,'fro') || ...
        abs(lam.I1)>1e-10*sqrt(lam.I0*lam.I2)
    error('Navier reference requires symmetric specially orthotropic laminate/inertia.');
end
al=m*pi/a; be=n*pi/b; D=lam.D; S=lam.S;
ref.K=[S(1,1)*al^2+S(2,2)*be^2 S(1,1)*al S(2,2)*be; ...
    S(1,1)*al D(1,1)*al^2+D(3,3)*be^2+S(1,1) (D(1,2)+D(3,3))*al*be; ...
    S(2,2)*be (D(1,2)+D(3,3))*al*be D(2,2)*be^2+D(3,3)*al^2+S(2,2)];
ref.M=diag([lam.I0 lam.I2 lam.I2]);
[V,L]=eig(ref.K,ref.M); [val,idx]=sort(real(diag(L))); V=V(:,idx);
for j=1:3, V(:,j)=V(:,j)/sqrt(V(:,j)'*ref.M*V(:,j)); end
ref.omega=sqrt(val); ref.frequencyHz=ref.omega/(2*pi); ref.vectors=V;
ref.force=[0;0;0];
if mod(m,2)==1 && mod(n,2)==1, ref.force(1)=16/(m*n*pi^2); end
end
