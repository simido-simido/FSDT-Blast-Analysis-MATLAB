function p = blast_pressure(t,pulse)
if pulse.pmax<=0 || pulse.td<=0 || pulse.decay<=0 || ...
        pulse.pmin<0 || pulse.tm<=0 || ...
        any(~isfinite([pulse.pmax pulse.td pulse.decay pulse.pmin pulse.tm]))
    error('Invalid pulse: positive pmax, td, decay, tm and nonnegative pmin required.');
end
p=zeros(size(t)); td=pulse.td;
positive=t>=0 & t<td;
switch lower(pulse.type)
    case {'friedlander','friedlander_cubic','extended'}
        x=t(positive)/td;
        p(positive)=pulse.pmax*(1-x).*exp(-pulse.decay*x);
    case 'triangular'
        p(positive)=pulse.pmax*(1-t(positive)/td);
    case 'rectangular'
        p(positive)=pulse.pmax;
    otherwise, error('Unknown pulse type %s.',pulse.type);
end
if pulse.includeNegative
    switch lower(pulse.type)
        case 'friedlander_cubic'
            neg=t>=td & t<td+pulse.tm; s=(t(neg)-td)/pulse.tm;
            p(neg)=-6.75*pulse.pmin*s.*(1-s).^2;
        case 'extended'
            neg=t>=td; x=t(neg)/td;
            p(neg)=pulse.pmax*(1-x).*exp(-pulse.decay*x);
    end
end
end
