function [A] = rand_phasor(nx,ny,h,nvp)
    % RAND_PHASOR Generate a random PhasorArray with a specified symmetry.
    %
    %   A = RAND_PHASOR(NX, NY, H, <name-value arguments>) creates a 3D array
    %   of complex phasors representing a time-dependent matrix with harmonic components.
    %
    %   Inputs:
    %     NX  - (integer) Number of rows.
    %     NY  - (integer) Number of columns.
    %     H   - (integer) Number of harmonics.
    %
    %   Name-Value Arguments:
    %     'symmetry' (string array) - Symmetry class of the result (default: "real").
    %                               Pass [] for a full complex array. The 14 names,
    %                               in signed pairs, are:
    %                                 symmetric      skewSymmetric
    %                                 real           imaginary
    %                                 even           odd
    %                                 hermitian      skewHermitian
    %                                 paraSymmetric  skewParaSymmetric
    %                                 paraConjugate  skewParaConjugate
    %                                 paraHermitian  skewParaHermitian
    %                               PHASORSYMMETRY states what each one means and
    %                               which combinations are reachable.
    %     'T' (scalar)            - Period of the system (default: 2*pi).
    %     'output' (char)         - Output format: 'NDarray' (default) or 'PhasorArray'.
    %     'average_power_decay'   - Decay rate of harmonic amplitudes (default: 2).
    %
    %   Outputs:
    %     A - (3D array or PhasorArray) Randomly generated PhasorArray.
    %
    %   Notes:
    %     - Generates phasors with complex Gaussian components, then projects onto
    %       the requested symmetry class.
    %     - The decay factor controls power reduction of higher-order harmonics.
    %     - 'symmetry' constrains structure, never the spectrum. For prescribed
    %       Floquet exponents use PhasorArray.randomWithNPole; for positive
    %       definiteness use PhasorArray.randomSPD.
    %     - A sparsity pattern is not a symmetry and has no name here. Use TRIU,
    %       TRIL or a product with a constant mask, none of which change the
    %       harmonic order and all of which preserve realness:
    %           U = triu(PhasorArray.random(3, 3, 5));       % A(t) upper triangular
    %
    %   Example
    %       A = rand_phasor(3, 3, 5, "symmetry", ["real" "symmetric"]);
    %       B = rand_phasor(3, 3, 5, "symmetry", []);          % full complex
    %
    %   See also: phasorSymmetry, PhasorArray, PhasorArray.randomSPD,
    %             PhasorArray.randomWithNPole, randomPhasorArrayWithPole.
arguments
    nx
    ny
    h
    nvp.symmetry string = "real"
    nvp.T=2*pi
    nvp.output {mustBeMember(nvp.output,{'NDarray','PhasorArray'})} ='NDarray' %or PhasorArray
    nvp.average_power_decay=2;
end

decay = nvp.average_power_decay;
A = pvalue(phasorSymmetry(PhasorArray(drawDecay(nx,ny,h,decay,nvp.T)), nvp.symmetry));

switch nvp.output
    case 'PhasorArray'
        A=PhasorArray(A);
    otherwise
end

end

%% =========================================================================
function A = drawDecay(nx,ny,h,decay,T)
% Unconstrained draw, harmonic amplitude decaying as 1/(|k|+1)^decay.
A = zeros(nx,ny,2*h+1);
for hh = 0:h
    w = 1/(hh+1)^decay/T;
    A(:,:,h+1+hh) = draw(nx,ny)*w;
    if hh > 0
        A(:,:,h+1-hh) = draw(nx,ny)*w;
    end
end
end

function M = draw(nx,ny)
M = ((rand(nx,ny)-0.5)*2 + 1i*(rand(nx,ny)-0.5)*2).*(rand(nx,ny)*8);
end
