function [fixed,free] = fsdt_boundary(mesh,bc)
%FSDT_BOUNDARY Explicit transverse and membrane constraints.
% ssss matches Ferreira: w=0 all edges, theta_x=0 on y edges,
% theta_y=0 on x edges. cccc: w=theta_x=theta_y=0 all edges.
n=mesh.nn; x=mesh.nodes(:,1); y=mesh.nodes(:,2);
tol=1e-10*max(mesh.a,mesh.b);
x0=find(abs(x)<tol); x1=find(abs(x-mesh.a)<tol);
y0=find(abs(y)<tol); y1=find(abs(y-mesh.b)<tol);
edge=unique([x0;x1;y0;y1]);
switch lower(bc.transverse)
    case 'ssss'
        tx=unique([y0;y1]); ty=unique([x0;x1]);
    case 'cccc'
        tx=edge; ty=edge;
    otherwise, error('transverse support must be ssss or cccc.');
end
switch lower(bc.membrane)
    case 'ferreira'
        if strcmpi(bc.transverse,'ssss'), u=x0; v=y0;
        else, u=edge; v=edge; end
    case 'immovable'
        u=edge; v=edge;
    case 'minimal'
        % Remove two in-plane translations and one rigid in-plane rotation.
        u=1; v=[1; mesh.nx+1];
    otherwise, error('membrane restraint must be ferreira, immovable, or minimal.');
end
fixed=unique([edge;tx+n;ty+2*n;u+3*n;v+4*n]);
free=setdiff((1:mesh.ndof)',fixed);
end
