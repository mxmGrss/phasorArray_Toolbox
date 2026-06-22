function out = ZPNSequence(dephase,order,include0)
    % ZPNSEQUENCE Produce PhasorArray representation of zero, positive, and negative sequence transformation.
    %  out = ZPNSEQUENCE(dephase, order) returns a PhasorArray representation of the zero, positive, and negative sequence transformation matrix.
    %
    %  Defining a = cos(θ) + j sin(θ), the zero, positive, and negative sequence transformation is defined as:
    %   P(θ) = [ 1 1 1;
    %            1 a^2 a;
    %            1 a a^2]
    %   where θ = k ϑ + dephase, with k the order of the transformation,
    %    and ϑ the angle of the transformation, often the electrical angle.
    %
    %  The zero, positive, and negative sequence transformation is used to transform three-phase quantities from abc to the ZPN frame.
    arguments
        dephase = 0
        order {mustBePositive, mustBeInteger} = 1
        include0=true
    end
    a = exp(2/3*1i*pi*order+dephase*1i);
    out=PhasorArray(cat(3,zeros(3,3,2),[1 1 1;1 a^2 a;1 a a^2]));
    if ~include0
        out = out{2:end,:}; %remove the zero sequence
    end
end
