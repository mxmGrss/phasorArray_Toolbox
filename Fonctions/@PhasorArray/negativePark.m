function out = negativePark(dephase,order,include0)
    % NEGATIVEPARK Produce PhasorArray representation of negative Park transformation.
    %  out = NEGATIVEPARK() returns a PhasorArray representation of the negative Park transformation matrix.
    %
    %  The negative dq0 transform is defined as P(θ) = negativeRot(θ) * Clark(),
    %    with θ = k ϑ + dephase. where k is the order of the transformation,
    %    and ϑ is the angle of the transformation, often the electrical angle.
    %
    %   Hence P(θ) = [cos(θ) cos(θ+2π/3) cos(θ-2π/3);
    %                 -sin(θ) -sin(θ+2π/3) -sin(θ-2π/3);
    %                 1/2     1/2         1/2      ]*2/3
    %
    %  The negative Park transformation is used to transform three-phase quantities from abc to the negative dq0 reference frame with park transform, conserving amplitude but not power.
    %
    %  See also: PARK, CLARK
    arguments
        dephase = 0
        order {mustBePositive, mustBeInteger} = 1
        include0=true
    end
    nRotPA = PhasorArray.negativeRotdq0(dephase, order);
    clarkPA = PhasorArray.Clark();
    out = nRotPA*clarkPA;
    if ~include0
        out = out{1:2,:}; %remove the zero sequence
    end
end

