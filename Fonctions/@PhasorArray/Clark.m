function out = Clark(include0)
    % CLARK Produce PhasorArray representation of Clark transformation.
    %  out = CLARK() returns a PhasorArray representation of the Clark transformation matrix.
    %
    %  The Clark transformation is defined as the matrix:
    %  Clark = (2/3) * [  1   -1/2  -1/2;
    %                     0   sqrt(3)/2  -sqrt(3)/2;
    %                    1/2   1/2   1/2 ]
    %
    %  The Clark transformation is used to transform three-phase quantities from abc to the αβ0 reference frame.
    %
    %  See also: PARK, DQ0, CONCORDIA

    arguments
        include0 = true
    end
    out = PhasorArray([1, -1/2, -1/2;
        0, sqrt(3)/2, -sqrt(3)/2;
        1/2, 1/2, 1/2]) * (2/3);
    if ~include0
        out = out{1:2,:}; %remove the zero sequence
    end
end


