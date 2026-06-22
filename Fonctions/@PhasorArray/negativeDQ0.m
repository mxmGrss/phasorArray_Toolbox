function out = negativeDQ0(dephase, order, include0)
    % NEGATIVEDQ0 Produce PhasorArray representation of negative dq0 transformation.
    %  out = NEGATIVEDQ0(dephase, order) returns a PhasorArray representation of the negative dq0 transformation matrix.
    %
    %  The negative dq0 transform is defined as P(θ) = negativeRotdq0(θ) * Concordia(),
    %    with θ = k ϑ + dephase, where k is the order of the transformation,
    %    and ϑ is the angle of the transformation, often the electrical angle.
    %
    %   Hence P(θ) = [cos(θ) cos(θ+2π/3) cos(θ-2π/3);
    %                 sin(θ) sin(θ+2π/3) sin(θ-2π/3);
    %                  1/√2     1/√2         1/√2      ]*√(2/3)
    %
    %  The negative dq0 transformation is used to transform three-phase quantities from abc to the negative dq0 reference frame.
    %
    %  See also: PARK, CLARK, CONCORDIA
    arguments
        dephase = 0
        order {mustBePositive, mustBeInteger} = 1
        include0 = true
    end
    nRotPA = PhasorArray.negativeRotdq0(dephase, order);
    concordiaPA = PhasorArray.Concordia();
    out = nRotPA * concordiaPA;
    if ~include0
        out = out{1:2, :}; % Remove the zero sequence
    end
end

