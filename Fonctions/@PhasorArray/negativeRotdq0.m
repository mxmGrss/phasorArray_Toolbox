function out = negativeRotdq0(dephase,order,include0)
    % NEGATIVEROTDQ0 Produce PhasorArray representation of negative rotation transformation.
    %  out = NEGATIVEROTDQ0(dephase, order) returns a PhasorArray representation of the rotation Park transformation matrix.
    %
    %  The negative rotation transformation is used to transform αβ0 quantities to the negative (indirect) dq0 reference frame.
    %
    %  It is defined as the matrix C(θ) = [ cos(θ)  -sin(θ) 0;
    %                                      -sin(θ)  -cos(θ) 0;
    %                                       0       0      1 ]
    %                   where θ = k ϑ + dephase, with k the order of the transformation,
    %                    and ϑ the angle of the transformation, often the electrical angle.
    %
    %   If include0 is false, the zero sequence is removed from the transformation, resulting in the 2x2 rotation matrix:
    %   R(θ) = [cos(θ)  sin(θ);
    %           sin(θ) -cos(θ)]
    %
    %  Inputs:
    %     dephase - (scalar, optional) Phase shift of the transformation. Default is 0.
    %     order   - (integer, optional) Harmonic order of the transformation. Default is 1.
    %     include0 - (logical, optional) Include zero sequence in the transformation. Default is true.
    %
    %  See also: PARK, DQ0
    arguments
        dephase = 0
        order {mustBePositive, mustBeInteger} = 1
        include0=true
    end
    cos = PhasorArray.cos(dephase, order);
    sin = PhasorArray.sin(dephase, order);
    out = [cos  -sin 0
        -sin -cos 0
        0 0 1];
    if ~include0
        out = out{1:2,1:2}; %remove the zero sequence
    end
end

