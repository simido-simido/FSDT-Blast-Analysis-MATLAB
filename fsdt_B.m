function [Bm,Bb,Bs] = fsdt_B(N,dN)
%FSDT_B Ferreira's Q4 strain-displacement matrices, preserved field ordering.
Bm=zeros(3,20); Bb=zeros(3,20); Bs=zeros(2,20);
Bb(1,5:8)=dN(:,1)'; Bb(2,9:12)=dN(:,2)';
Bb(3,5:8)=dN(:,2)'; Bb(3,9:12)=dN(:,1)';
Bm(1,13:16)=dN(:,1)'; Bm(2,17:20)=dN(:,2)';
Bm(3,13:16)=dN(:,2)'; Bm(3,17:20)=dN(:,1)';
Bs(1,1:4)=dN(:,1)'; Bs(2,1:4)=dN(:,2)';
Bs(1,5:8)=N'; Bs(2,9:12)=N';
end
