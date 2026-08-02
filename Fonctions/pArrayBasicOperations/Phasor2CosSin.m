%{
Phasor2CosSin - Converts a phasor representation to cosine-sine representation.

This function takes a matrix or higher-dimensional array of phasors and converts 
it into its corresponding cosine-sine representation. The function supports 
various options for handling the input data, including specifying the harmonic 
storage dimension, restricting to positive harmonics, and handling pre-converted 
cosine-sine data.

Syntax:
    [Mcs, hM] = Phasor2CosSin(Mph, nvp)

Inputs:
    Mph - Input matrix or array representing the phasors. The dimensions and 
          interpretation depend on the optional arguments.

    nvp - optional name-value pairs:
        only_posk (logical, default: false) - If true, restricts the output to 
            positive harmonics only.
        harmodim (numeric, default: []) - Specifies the dimension along which 
            harmonics are stored. If empty, a warning is issued for 2D inputs.
        already_cs (logical, default: false) - If true, assumes the input is 
            already in cosine-sine format and skips conversion.

Outputs:
    Mcs - Output matrix or array in cosine-sine representation. The dimensions 
          depend on the input and options provided.
    hM  - The number of harmonics (half the size of the third dimension minus one).

Warnings:
    - If a 2D input is provided without specifying 'harmodim', the function 
      assumes it represents a 0th phasor of a constant matrix and issues a warning.

Notes:
    - The function supports multi-dimensional arrays, with the 4th and 5th 
      dimensions used for additional indexing.
    - The cosine-sine representation is constructed by combining the real and 
      imaginary parts of the phasors.

Example:
    % Example usage of Phasor2CosSin
    Mph = rand(3, 3, 5); % Random phasor matrix
    [Mcs, hM] = Phasor2CosSin(Mph, struct('only_posk', true));

%}
function [Mcs,hM] = Phasor2CosSin(Mph,nvp)
arguments
    Mph
    nvp.only_posk=false
    nvp.harmodim=[]
    nvp.already_cs=false
end

if ismatrix(Mph)
    if isempty(nvp.harmodim)
        warning('Phasor2CosSin:2DInputAmbiguous', '2D input interpreted as the 0th phasor of a constant matrix. Specify harmodim if the input is a vector of harmonics.')
    elseif nvp.harmodim==2
        Mph=permute(Mph,[1 3 2]);
    elseif nvp.harmodim==1
        Mph=permute(Mph,[2 3 1]);
    end
end

if nvp.already_cs
    Mcs = Mph;
    hM = (size(Mph,3)-1)/2;
    return
end

n4=size(Mph,4);
n5=size(Mph,5);
for iter_ii=1:n5
    for iter_i=1:n4
        if ~nvp.only_posk
            Mph_t=Mph(:,:,(end+1)/2:end,iter_i,iter_ii);
        end
        hM=size(Mph_t,3)-1;
        Mcs(:,:,:,iter_i,iter_ii)=cat(3,-2*imag(flip(Mph_t(:,:,2:end),3)), real(Mph_t(:,:,1)), 2*real(Mph_t(:,:,2:end)));
    end
end
end

