function [Mt,t] = PhasorArray2time(Mph,T,t,nvp)
% PHASORARRAY2TIME Evaluate a periodic matrix function A(t) or A(θ) from its phasor representation.
%
%   PHASORARRAY2TIME reconstructs the time-domain representation of a periodic 
%   matrix function A(t) given its phasor coefficients. It can evaluate:
%   - A(t) for a fixed period `T`
%   - A(θ) for a vector of angles `θ`
%   - A(t) when `T` and `t` are both vectors, in which case `T` is interpreted as θ(t).
%
%   Inputs:
%     Mph  - (m × n × (2h+1) array) Phasor coefficients of the periodic matrix A(t).
%     T    - (double or vector, optional) Period of the function or a vector of angles `θ`:
%              - If `T` is a **scalar**, it represents the period of A(t).
%              - If `T` is a **vector**:
%                 - If `t` is empty, it is interpreted as `θ`, and A(θ) is computed.
%                 - If `T` and `t` have matching sizes, `T` is interpreted as a time-dependent frequency phase `θ(t)`.
%              - Default: `2*pi`.
%     t    - (vector, optional) Time instants at which to evaluate A(t).
%              - If `T` is a vector and `t` is empty, then `T` is interpreted as `θ`, and `t = T`.
%              - If `T` and `t` have matching sizes, `T` is interpreted as `θ(t)`, defining a time-varying phase.
%              - If empty, a default time grid is computed.
%     nvp  -  Name-value pair arguments:
%              - 'plot' (logical): Plot matrix coefficients versus time (default: `false`).
%              - 'explosed' (logical): Plot each coefficient in a separate subplot (default: `true`).
%              - 'hold' (logical): Hold the current plot (default: `false`).
%              - 'DispImag' (logical): Display the imaginary part of the matrix (default: `false`).
%              - 'DispReal' (logical): Display the real part of the matrix (default: `true`).
%              - 'ZeroCentered' (logical): Center the plot around zero (default: `false`).
%              - 'forceReal' (logical): Force the matrix to be real-valued (default: `[]`).
%              - 'checkReal' (logical): Check if the matrix is real-valued (default: `false`).
%              - 'checkRealTol' (double): Tolerance for real-valued check (default: `1e-8`).
%              - 'title' (string): Title of the plot (default: `[]`).
%              - 'plot3D' (logical): Plot the matrix coefficients in 3D (default: `false`).
%              - 'LineStyle' (string): Line type for the plot (default: `'-'`).
%              - 'GlobalYLim' (logical): Set uniform y-limits for all subplots (default: `false`).
%              - 'linkaxes' (string): Link axes ('x', 'y', 'xy', etc.) (default: `'x'`).
%              - 'computationMethod' (integer): Method used for evaluation (default: `1`).
%              - 'providedPhasorForm' (string): Phasor representation ('exp' or 'SinCos') (default: `"exp"`).
%              - 'parent' (handle): Parent figure or axes (default: `[]`).
%
%   Outputs:
%     Mt   - (m × n × length(t) array) Evaluated matrix A(t) or A(θ).
%     t    - (vector) Time instants (or angles θ) at which A(t) is evaluated.
%
%   Behavior:
%     - **Evaluation Modes**:
%         - If `T` is a scalar and `t` is provided, computes **A(t) for a fixed period T**.
%         - If `T` is a vector and `t` is empty, computes **A(θ)** using `θ = T`.
%         - If `T` and `t` are both vectors of the same length, computes **A(t) with θ(t) as a time-varying frequency phase**.
%     - **Automatic Time Grid**:
%         - If `t` is empty, it is computed automatically based on `T` and the highest harmonic.
%     - **Real-Valued Constraints**:
%         - Allows optional enforcement of real-valued output.
%     - **Plotting Options**:
%         - Can generate time-domain plots of individual matrix elements.
%
%   Example Usage:
%     % Evaluate a periodic matrix at given time points
%     Mph = rand(3, 3, 5);
%     T = 2*pi;
%     t = 0:0.01:1;
%     [Mt, t] = PhasorArray2time(Mph, T, t);
%
%     % Evaluate at a vector of angles instead of time
%     theta = linspace(0, 2*pi, 100);
%     [Mt, theta] = PhasorArray2time(Mph, theta);
%
%     % Evaluate A(t) with a time-varying phase θ(t)
%     t = linspace(0, 1, 100);
%     theta_t = 2*pi*t.^2;  % Example of a quadratic phase modulation
%     [Mt, t] = PhasorArray2time(Mph, theta_t, t);
%
%     % Generate a time-domain plot
%     [Mt, t] = PhasorArray2time(Mph, T, t, 'plot', true);
%
%   See also: timeArray2Phasors, stemPhasor.

arguments
    Mph
    T=2*pi
    t=[]
    nvp.plot logical         = false
    nvp.explosed logical     = true
    nvp.hold logical         = false
    nvp.DispImag logical     = false
    nvp.DispReal logical     = true
    nvp.ZeroCentered logical = false
    nvp.forceReal logical    = []
    nvp.checkReal logical    = false
    nvp.checkRealTol         = 1e-8
    nvp.title                = []
    nvp.plot3D logical       = false
    nvp.LineStyle             = '-'
    nvp.GlobalYLim logical   = false
    nvp.linkaxes             = 'x';
    nvp.computationMethod    = 1
    nvp.providedPhasorForm {mustBeMember(nvp.providedPhasorForm,["exp","SinCos"])} = "exp"
    nvp.parent               = [];
    nvp.grid {mustBeMember(nvp.grid,["off","on","minor"])}                = 'on'
    nvp.squeeze logical      = false
end

% No release gate here any more: methods 1 and 2 went through tensorprod, which
% is R2022a, so older releases were pushed onto the slower sparse method 3. Both
% now go through harmonicCombine and a plain product, so every release takes the
% fast path.



if isempty(nvp.forceReal)
    if nvp.providedPhasorForm == "exp"
        if isrealp(Mph)
            nvp.forceReal = true;
        end
    else
        nvp.forceReal = false;
    end
end

% User specified that the input PhasorArray is in SinCos form
if nvp.providedPhasorForm == "SinCos"
    % Ensure that the matrix is real valued
    if ~nvp.forceReal
        nvp.forceReal = true;
        warning('PhasorArray2time:sinCosForceReal', 'SinCos phasor form is only compatible with real-valued matrices. Forcing real-valued computation.')
    end

    % Convert the matrix to SinCos form if it is not already
    if isa(Mph, 'PhasorArray')
        Mph = Mph.SinCosForm();
        warning('PhasorArray2time:sinCosConversion', 'Input PhasorArray converted to SinCos form.')
    end


end

    % Get dimensions
    if isa(Mph, 'PhasorArray')
        Mph_val = Mph.value;
    else
        Mph_val = Mph;
    end
    
    nx = size(Mph_val, 1);
    ny = size(Mph_val, 2);
    h_len = size(Mph_val, 3);
    h = (h_len - 1) / 2;

% An empty t means T carries the phase samples, so the abscissa is an angle.
% Decided here because t is resolved just below and the distinction is lost.
if isempty(t)
    nvp.xlabelStr = 'angle (rad)';
else
    nvp.xlabelStr = 'time (sec)';
end

% Compute the time vector if not provided
if and(numel(t) < 3 , numel(T) == 1)
    dt = computeTimeStep(T, h);
    t = computeTimeVector(t, dt, T);
end

% Ensure T and t are row vectors
T = reshape(T, 1, []);
t = reshape(t, 1, []);

% Compute theta
theta = computeTheta(T, t);

% Compute the basis for evaluation
if ~nvp.forceReal
    % Complex-valued basis
    eit = exp(1i * (-h:h)' * theta);
    Meval = Mph_val;
else
    % Real-valued basis (Sin/Cos representation)
    if nvp.providedPhasorForm == "exp"
        % Convert exponential phasors to real conjugate-symmetric form
        Mphr  = real(Mph_val + flip(Mph_val, 3)) / 2 + 1i * imag(Mph_val - flip(Mph_val, 3)) / 2;
        Meval = real(cat(3, 1i * (flip(Mphr(:, :, h+2:end), 3) - Mphr(:, :, 1:h)), Mphr(:, :, h+1), (Mphr(:, :, h+2:end) + flip(Mphr(:, :, 1:h), 3))));
    else
        Meval = Mph_val;
    end
    eit = [sin((h:-1:1)' * theta); cos((0:h)' * theta)];
end

if isa(t,'sym')
    nvp.computationMethod=3;
end

%choose the computation method to compute the matrix at each time
switch nvp.computationMethod
    case 1
        %best time but tensorprod reliant so matlab >22a
        Mt=harmonicCombine(Meval,double(eit)); %est un 3D array dont Mt(:,:,k) est M(t(k))

    case 2
        %bad time very slow
        reM=reshape(Meval,nx,[],1);
        reEit=kron(eit,eye(ny));
        rMt=reM*reEit;   % 2-D operands: the mode-2 contraction is a plain product
        Mt=reshape(rMt,nx,ny,[]);

    otherwise
        %a bit better thanks to sparse
        reM=reshape(Meval,nx,[],1);
        reEit=kron(eit,speye(ny));
        rbMt=reM*reEit;
        Mt=reshape(rbMt,nx,ny,[]);
end

% Check if the imaginary part of Mt is negligible
if nvp.checkReal
    if all(abs(imag(Mt)) < nvp.checkRealTol, 'all')
        Mt = real(Mt);
    else
        warning('PhasorArray2time:notReal', 'Output is complex-valued despite checkReal flag (imaginary part above tolerance).');
    end
end


% Plotting lives in TimeArray2plot: most callers here only want [Mt, t].
if nvp.plot
    t = TimeArray2plot(Mt, t, T, nvp);
end

if nvp.squeeze
    Mt = squeeze(Mt);
end

    function Mph = convertToNumeric(Mph)
        if isa(Mph, 'PhasorArray')
            Mph = Mph.Value;
        elseif isa(Mph, "ndsdpvar") || isa(Mph, 'sdpvar')
            Mph = value(Mph);
        elseif isa(Mph, 'sym')
            Mph = vpa(value(Mph));
            nvp.computationMethod = 3;
        end
    end
    function dt = computeTimeStep(T, h)
        % This function computes the time step based on the period T and the harmonic order h
        % Inputs: T - period, h - harmonic order
        % Output: dt - time step

        dt = T / 2 ^ (nextpow2(max((h * 8), 64)));
    end

    function t = computeTimeVector(t, dt, T)
        % This function computes the time vector based on the input t, the time step dt, and the period T
        % Inputs: t - input time, dt - time step, T - period
        % Output: t - time vector

        switch numel(t)
            case 0 %default one period evaluation
                t = 0:dt:T-dt;
            case 1
                % DO NOTHING scalar case of evaluation at only one given time
            case 2
                t = t(1):dt:t(2);
        end
    end

    function theta = computeTheta(T, t)
        % This function computes the theta value based on the period T and the time vector t
        % Inputs: T - period, t - time vector
        % Output: theta - computed theta value

        if numel(T) > 1
            theta = T;
        else
            theta = 2 * pi / T * t;
        end
    end
end
