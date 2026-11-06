function tr = fsdt_newmark(M,C,K,funit,t,p,u0,v0)

t=t(:)'; p=p(:)'; nf=size(M,1); nt=numel(t);
if nt<2 || numel(p)~=nt || any(~isfinite(t)) || any(~isfinite(p))
    error('Finite t,p with equal lengths and at least two samples required.');
end
dt=t(2)-t(1);
if dt<=0 || max(abs(diff(t)-dt))>1e-8*dt
    error('Newmark requires a positive uniform time step.');
end
if nargin<8 || isempty(u0), u0=zeros(nf,1); end
if nargin<9 || isempty(v0), v0=zeros(nf,1); end
if numel(u0)~=nf || numel(v0)~=nf, error('Initial-state dimensions disagree.'); end
beta=1/4; gamma=1/2;
a0=1/(beta*dt^2); a1=gamma/(beta*dt); a2=1/(beta*dt);
a3=1/(2*beta)-1; a4=gamma/beta-1; a5=dt*(gamma/(2*beta)-1);
U=zeros(nf,nt); V=U; A=U; U(:,1)=u0(:); V(:,1)=v0(:);
A(:,1)=M\(funit*p(1)-C*V(:,1)-K*U(:,1));
Keff=K+a0*M+a1*C;
[R,flag]=chol(Keff);
if flag~=0, error('Effective stiffness is not positive definite.'); end
Ek=zeros(1,nt); Es=Ek; work=Ek; diss=Ek; residual=Ek;
Ek(1)=0.5*V(:,1)'*M*V(:,1); Es(1)=0.5*U(:,1)'*K*U(:,1);
for j=2:nt
    rhs=funit*p(j)+M*(a0*U(:,j-1)+a2*V(:,j-1)+a3*A(:,j-1)) ...
        +C*(a1*U(:,j-1)+a4*V(:,j-1)+a5*A(:,j-1));
    U(:,j)=R\(R'\rhs);
    A(:,j)=a0*(U(:,j)-U(:,j-1))-a2*V(:,j-1)-a3*A(:,j-1);
    V(:,j)=V(:,j-1)+dt*((1-gamma)*A(:,j-1)+gamma*A(:,j));
    Ek(j)=0.5*V(:,j)'*M*V(:,j); Es(j)=0.5*U(:,j)'*K*U(:,j);
    du=U(:,j)-U(:,j-1); vavg=0.5*(V(:,j)+V(:,j-1));
    work(j)=work(j-1)+0.5*(p(j)+p(j-1))*(funit'*du);
    diss(j)=diss(j-1)+dt*(vavg'*C*vavg);
    r=M*A(:,j)+C*V(:,j)+K*U(:,j)-funit*p(j);
    scale=norm(M*A(:,j))+norm(C*V(:,j))+norm(K*U(:,j))+norm(funit*p(j));
    residual(j)=norm(r)/max(scale,realmin);
end
tr.time=t; tr.pressure=p; tr.U=U; tr.V=V; tr.A=A;
tr.kinetic=Ek; tr.strain=Es; tr.work=work; tr.dissipated=diss;
tr.energyError=Ek+Es+diss-work-Ek(1)-Es(1);
tr.equilibriumResidual=residual; tr.dt=dt;
end
