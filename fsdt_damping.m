function [C,coeff] = fsdt_damping(model,modal,damping)
%FSDT_DAMPING Rayleigh damping fitted at two explicitly selected flexural modes.
z=damping.ratios(:);
if numel(z)==1, z=[z;z]; end
if numel(z)~=2 || any(z<0) || any(~isfinite(z))
    error('damping.ratios must contain one or two nonnegative ratios.');
end
if all(z==0), C=sparse(size(model.Mff,1),size(model.Mff,2)); coeff=[0;0]; return; end
idx=damping.modes;
if numel(idx)~=2 || any(idx<1) || any(idx~=round(idx)) || ...
        max(idx)>numel(modal.flexural) || idx(1)==idx(2)
    error('Choose two available distinct flexural mode indices for damping.');
end
w=modal.omega(modal.flexural(idx));
if abs(w(1)-w(2))<1e-8*max(w), error('Damping fit needs distinct frequencies.'); end
coeff=[1./(2*w), w/2]\z;
if any(coeff<0), error('Requested ratios produce negative Rayleigh coefficients.'); end
C=coeff(1)*model.Mff+coeff(2)*model.Kff;
end
