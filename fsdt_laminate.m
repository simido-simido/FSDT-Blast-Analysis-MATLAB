function lam = fsdt_laminate(plies, shearCorrection)
if isempty(plies), error('At least one ply is required.'); end
kap = shearCorrection(:);
if numel(kap)==1, kap=[kap;kap]; end
if numel(kap)~=2 || any(~isfinite(kap)) || any(kap<=0)
    error('shearCorrection must be positive scalar or [k_xz k_yz].');
end
n = numel(plies);
fields = {'t','E1','E2','G12','G13','G23','rho'};
for k=1:n
    for j=1:numel(fields)
        v=plies(k).(fields{j});
        if ~isscalar(v) || ~isfinite(v) || v<=0
            error('Invalid positive ply property %s in ply %d.',fields{j},k);
        end
    end
    if ~isscalar(plies(k).angle) || ~isfinite(plies(k).angle) || ...
            ~isscalar(plies(k).nu12) || ~isfinite(plies(k).nu12)
        error('Invalid angle/nu12 in ply %d.',k);
    end
end
lam.h=sum([plies.t]);
lam.z=[-lam.h/2 -lam.h/2+cumsum([plies.t])];
lam.A=zeros(3); lam.B=zeros(3); lam.D=zeros(3); lam.S=zeros(2);
lam.Qbar=zeros(3,3,n); lam.Qs=zeros(2,2,n);
lam.I0=0; lam.I1=0; lam.I2=0; lam.plies=plies;
for k=1:n
    p=plies(k); nu21=p.nu12*p.E2/p.E1; den=1-p.nu12*nu21;
    if den<=0, error('Non-positive orthotropic stiffness denominator.'); end
    q11=p.E1/den; q22=p.E2/den; q12=p.nu12*p.E2/den; q66=p.G12;
    c=cos(p.angle*pi/180); s=sin(p.angle*pi/180);
    q=zeros(3);
    q(1,1)=q11*c^4+2*(q12+2*q66)*c^2*s^2+q22*s^4;
    q(2,2)=q11*s^4+2*(q12+2*q66)*c^2*s^2+q22*c^4;
    q(1,2)=(q11+q22-4*q66)*c^2*s^2+q12*(c^4+s^4);
    q(1,3)=(q11-q12-2*q66)*c^3*s-(q22-q12-2*q66)*c*s^3;
    q(2,3)=(q11-q12-2*q66)*c*s^3-(q22-q12-2*q66)*c^3*s;
    q(3,3)=(q11+q22-2*q12-2*q66)*c^2*s^2+q66*(c^4+s^4);
    q(2,1)=q(1,2); q(3,1)=q(1,3); q(3,2)=q(2,3);
    qs=[p.G13*c^2+p.G23*s^2, (p.G13-p.G23)*c*s; ...
        (p.G13-p.G23)*c*s, p.G13*s^2+p.G23*c^2];
    z0=lam.z(k); z1=lam.z(k+1);
    dz=z1-z0; dz2=(z1^2-z0^2)/2; dz3=(z1^3-z0^3)/3;
    lam.Qbar(:,:,k)=q; lam.Qs(:,:,k)=qs;
    lam.A=lam.A+q*dz; lam.B=lam.B+q*dz2; lam.D=lam.D+q*dz3;
    R=diag(sqrt(kap)); lam.S=lam.S+R*qs*R*dz;
    lam.I0=lam.I0+p.rho*dz;
    lam.I1=lam.I1+p.rho*dz2;
    lam.I2=lam.I2+p.rho*dz3;
end
lam.inertia=diag([lam.I0 lam.I2 lam.I2 lam.I0 lam.I0]);
lam.inertia(2,4)=lam.I1; lam.inertia(4,2)=lam.I1;
lam.inertia(3,5)=lam.I1; lam.inertia(5,3)=lam.I1;
lam.shearCorrection=kap;
end
