function out = zeroPosNegSequenceDQ(dephase,order,include0)
    % ZEROPOSNEGSEQUENCEDQ Produce PhasorArray representation of zero, positive, and negative sequence dq transformation.
    %  out = ZEROPOSNEGSEQUENCEDQ(dephase, order) returns a PhasorArray representation of the zero, positive, and negative sequence dq transformation matrix.
    %
    %  The zero, positive, and negative sequence dq transform is defined as:
    %    P(θ) = [zero(θ) ; dq(θ) ; negativeDQ(θ)],
    %    with θ = k ϑ + dephase, where k is the order of the transformation,
    %    and ϑ is the angle of the transformation, often the electrical angle.
    %
    %   Hence P(θ) = [  1/√2     1/√2         1/√2  ;
    %                 cos(θ) cos(θ-2π/3) cos(θ+2π/3);
    %                 sin(θ) sin(θ-2π/3) sin(θ+2π/3);
    %                 cos(θ) cos(θ+2π/3) cos(θ-2π/3);
    %                 sin(θ) sin(θ+2π/3) sin(θ-2π/3)]*√(2/3)
    %
    %  The zero, positive, and negative sequence dq transformation is used to transform three-phase quantities from abc to the dq0 reference frame.
    %
    %  See also: PARK, CLARK, DQ0, NEGATIVEDQ0

    arguments
        dephase = 0
        order {mustBePositive, mustBeInteger} = 1
        include0=true
    end
    dqPos = PhasorArray.dq0(dephase, order,false);
    dqNeg = PhasorArray.negativeDQ0(dephase, order,true);
    out = [dqNeg{3,:} ; dqPos ; dqNeg{1:2,:}]; %zero, positive dq, negative dq
    if ~include0
        out = out{2:end,:}; %remove the zero sequence
    end
end

