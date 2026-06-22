function out = Park(dephase, order,include0)
    % PARK Produce PhasorArray representation of PARK-dq0 transformation.
    %  out = PARK() returns a PhasorArray representation of the PARK transformation matrix.
    %
    %  The dq0 transform is defined as P(θ) = Park(θ) * Clark(),
    %    with θ = k ϑ + dephase. where k is the order of the transformation,
    %    and ϑ is the angle of the transformation, often the electrical angle.
    %
    %   Hence P(θ) = [cos(θ) cos(θ-2π/3) cos(θ+2π/3);
    %                 sin(θ) sin(θ-2π/3) sin(θ+2π/3);
    %                  1/√2     1/√2         1/√2      ]*√(2/3)
    %
    %  The PARK transformation is used to transform three-phase quantities from abc to the PARK dq0 reference frame, ocnserving amplitude of signals (but not power).
    %
    %  See also: PARK, CLARK
    arguments
        dephase = 0
        order {mustBePositive, mustBeInteger} = 1
        include0=true
    end
    rotPA = PhasorArray.Rotdq0(dephase, order);
    clarkPA = PhasorArray.Clark();
    out = rotPA*clarkPA;
    if ~include0
        out = out{1:2,:}; %remove the zero sequence
    end
end


