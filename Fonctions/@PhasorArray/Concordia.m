function out = Concordia(include0)
    % CONCORDIA Produce PhasorArray representation of power invariant Concordia transformation.
    %  out = CONCORDIA() returns a PhasorArray representation of the Concordia transformation matrix.
    %
    %  The CONCORDIA transformation is defined as the matrix : Concordia =     [   1      -1/2     -1/2  ;
    %                                                                      0      √3/2    -√3/2  ;
    %                                                                   1/√2      1/√2     1/√2  ]* √(2/3)
    %
    %  The CONCORDIA transformation is used to transform three-phase quantities from abc to the αβ0 reference frame, where power is preserved
    %
    %  See also: PARK, DQ0, CLARK
    arguments
        include0=true
    end
    out = PhasorArray([1, -1/2, -1/2;
        0, sqrt(3)/2, -sqrt(3)/2
        1/sqrt(2), 1/sqrt(2), 1/sqrt(2)])*sqrt(2/3);
    if ~include0
        out = out{1:2,:}; %remove the zero sequence
    end
end

%the phasor array representing Clark transformation, that preserves amplitude of signals but not power
