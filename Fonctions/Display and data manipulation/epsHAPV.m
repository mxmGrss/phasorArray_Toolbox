function [epsilons,eepsilon,adbParam] = epsHAPV(th, om, firstRotIsNan)
switch nargin
    case 2
        % Handle case when no firstRotIsNan is provided
        firstRotIsNan = true;
    otherwise
        firstRotIsNan = false;
end
[k, ~, ~] = find2piAntecedant(th,om,firstRotIsNan);

epsilons = zeros(size(k));
eepsilon = zeros(size(k));

%adiabatic parameter dot omega/ omega^2 = d/dphi omega / omega
%perform numerical differentiation using gradient
dotOmega = gradient(om, th); % Numerical differentiation to get dot omega
dotOmegaOverOmegaSquared = dotOmega ./ (om ); % Compute dot omega / omega²
adbParam = abs(dotOmegaOverOmegaSquared * 2*pi);

for ki = 1:numel(k)
    if isnan(k(ki))
        epsilons(ki) = NaN;
        eepsilon(ki) = NaN;
        continue
    end
    indexPReviousRev = (k(ki)):ki ; %all index of previous rotation
    omt = om(ki); %omega courant
    omtau = om(indexPReviousRev) ; %omega(tau) sur la derniere rot
    epsilons(ki) = max(abs((omtau - omt)./omtau)); %epsilon(t)

    %epsilon over the last rotation
    epstau = epsilons(indexPReviousRev) ; %epsilon sur la derniere rot
    %supp over last rotation
    eepsilon(ki) = max(abs(epstau)); %max de epsilon sur une rotation

end
end