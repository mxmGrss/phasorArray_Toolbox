function [k, f, ushifted] = find2piAntecedant(th,u)
%FIND2PIANTECEDANT Find 2π antecedent indices and compute phase-shifted signals.
%
%   [K, F] = FIND2PIANTECEDANT(TH) finds for each element TH(i), the largest 
%   index K(i) such that TH(K(i)) <= TH(i) - 2π, along with an interpolation 
%   factor F(i) for linear interpolation between TH(K(i)) and TH(K(i)+1).
%
%   [K, F, USHIFTED] = FIND2PIANTECEDANT(TH, U) additionally computes 
%   USHIFTED(theta) = U(theta - 2π), i.e., the signal values that occurred 
%   2π earlier in phase. This is particularly useful for analyzing periodic 
%   signals and their phase relationships, delay effects, and harmonic content.
%
%   This function is essential for phase unwrapping, periodic signal analysis,
%   and computing phase-delayed versions of signals for control and filtering
%   applications.
%
%   Inputs:
%       TH - Strictly increasing vector of phase values (numeric vector)
%       U  - (Optional) Signal values U(TH) corresponding to TH phases 
%            (numeric vector, same length as TH). If provided, USHIFTED 
%            will compute U(TH - 2π) through interpolation.
%
%   Outputs:
%       K        - Index vector where K(i) is the largest index such that 
%                  TH(K(i)) <= TH(i) - 2π. If no such index exists, K(i) = 1.
%       F        - Interpolation factor vector where F(i) represents the fractional
%                  position between TH(K(i)) and TH(K(i)+1) for the value TH(i) - 2π.
%                  If K(i) = 1 and no valid antecedent exists, F(i) = 0.
%       USHIFTED - Signal values at phase TH - 2π, i.e., USHIFTED(i) = U(TH(i) - 2π).
%                  Computed via linear interpolation: U(K(i)) + F(i) * (U(K(i)+1) - U(K(i))).
%                  Empty array [] if U is not provided or empty.
%
%   Mathematical relationships:
%       Phase antecedent: TH(i) - 2π ≈ TH(K(i)) + F(i) * (TH(K(i)+1) - TH(K(i)))
%       Signal at shifted phase: USHIFTED(i) = U(TH(i) - 2π) ≈ U(K(i)) + F(i) * (U(K(i)+1) - U(K(i)))
%
%   Examples:
%       % Analyze harmonic relationship in a periodic signal
%       th = linspace(0, 4*π, 200);
%       u = sin(th) + 0.3*sin(3*th);  % Fundamental + 3rd harmonic
%       [k, f, u_delayed] = find2piAntecedant(th, u);
%       
%       % Compare current signal with 2π-delayed version
%       figure;
%       plot(th, u, 'b-', th, u_delayed, 'r--');
%       legend('U(θ)', 'U(θ-2π)', 'Location', 'best');
%       xlabel('Phase θ'); ylabel('Signal Value');
%       title('Signal and its 2π Phase-Delayed Version');
%       
%       % Useful for control applications - compare current and previous cycle
%       th_control = linspace(0, 6*π, 300);
%       control_signal = square(th_control) + 0.1*randn(size(th_control));
%       [~, ~, prev_cycle] = find2piAntecedant(th_control, control_signal);
%       error_signal = control_signal - prev_cycle;  % Cycle-to-cycle error
%
%   Applications:
%       - Periodic control systems (comparing current vs. previous cycle)
%       - Harmonic analysis (phase relationships between harmonics)
%       - Signal processing (delay-based filtering and analysis)
%       - Vibration analysis (comparing successive periods)
%
%   See also: UNWRAP, ANGLE, MOD, INTERP1, CIRCSHIFT
%
%   PhasorArray Toolbox
%   Author: Maxime Grosso
%   Version: 1.0

    % Input validation
    arguments
        th {mustBeNumeric, mustBeVector, mustBeIncreasing}
        u = []
    end
    
    % Ensure th is a row vector for consistent indexing
    th = th(:)';
    n = length(th);
    
    % Initialize output arrays
    k = zeros(size(th));
    f = zeros(size(th));
    
    % Early return for empty or single-element input
    if n <= 1
        if n == 1
            k(1) = 1;
            f(1) = 0;
        end
        return;
    end
    
    last_k = 1; % Start from the first index
    
    for ii = 1:n
        target_value = th(ii) - 2*pi;
        
        % Continue the search from the last found index
        % Find the largest index where th(index) <= target_value
        while last_k < n && th(last_k) <= target_value
            last_k = last_k + 1;
        end
        
        % Adjust last_k to point to the largest valid index
        last_k = max(last_k - 1, 1);
        
        if last_k > 0 && th(last_k) <= target_value
            k(ii) = last_k;
            % Compute the fractional part f(ii) for interpolation
            if last_k < n
                f(ii) = (target_value - th(last_k)) / (th(last_k+1) - th(last_k));
            else
                f(ii) = 0; % At the boundary, no interpolation needed
            end
        else
            % If no such k(ii) exists, set k(ii) and f(ii) to default values
            k(ii) = 1;
            f(ii) = 0;
        end
    end

    if nargin > 1 && ~isempty(u)
        u= u(:)'; % Ensure u is a row vector
        % If a second input u is provided, compute the shifted values
        ushifted = zeros(size(u));
        for ii = 1:n
            ushifted(ii) = u(k(ii)) + f(ii) * (u(k(ii)+1) - u(k(ii)));
        end
    else
        ushifted = [];
    end
end

function mustBeIncreasing(th)
    % Validate that th is strictly increasing
    if ~issorted(th, 'strictascend')
        error('find2piAntecedant:NotIncreasing', ...
              'Input vector TH must be strictly increasing.');
    end
end
