function [epsilons,eepsilon,adbParam] = epsHAPV(th, om)
[k, ~, ~] = find2piAntecedant(th,om);

epsilons = zeros(size(k));
eepsilon = zeros(size(k));

%adiabatic parameter dot omega/ omega^2 = d/dphi omega / omega
%perform numerical differentiation using gradient
dotOmega = gradient(om, th); % Numerical differentiation to get dot omega
dotOmegaOverOmegaSquared = dotOmega ./ (om ); % Compute dot omega / omega²
adbParam = abs(dotOmegaOverOmegaSquared * 2*pi);

for ki = 1:numel(k)
    indexPReviousRev = (k(ki)):ki ;
    omt = om(ki);
    omtau = om(indexPReviousRev) ;
    epsilons(ki) = max(abs((omtau - omt)./omtau));

    %epsilon over the last rotation
    epstau = epsilons(indexPReviousRev) ;
    %supp over last rotation
    eepsilon(ki) = max(abs(epstau));

end
end