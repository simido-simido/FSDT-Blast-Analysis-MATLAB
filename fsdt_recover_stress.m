function stress = fsdt_recover_stress(lam,probe,displacements)

[Bm,Bb,Bs]=fsdt_B(probe.N,probe.dN);
ue=displacements(probe.dof,:); e0=Bm*ue; curvature=Bb*ue;
nt=size(ue,2); np=numel(lam.plies);
stress.globalPa=zeros(3,2,np,nt); stress.localPa=stress.globalPa;
stress.epsilon0=e0; stress.curvature=curvature; stress.shearStrain=Bs*ue;
for k=1:np
    c=cos(lam.plies(k).angle*pi/180); s=sin(lam.plies(k).angle*pi/180);
    T=[c*c s*s 2*c*s; s*s c*c -2*c*s; -c*s c*s c*c-s*s];
    for face=1:2
        z=lam.z(k+face-1);
        sig=lam.Qbar(:,:,k)*(e0+z*curvature);
        stress.globalPa(:,face,k,:)=reshape(sig,[3 1 1 nt]);
        stress.localPa(:,face,k,:)=reshape(T*sig,[3 1 1 nt]);
    end
end
stress.point=probe.point; stress.element=probe.element;
end
