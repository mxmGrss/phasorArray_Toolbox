function obj= ScalarPhasorArray(varargin,varg)
    % ScalarPhasorArray Creates a scalar phasor array from input data.
    %
    %   obj = ScalarPhasorArray(varargin, varg) creates a scalar PhasorArray, ie a scalar periodic signal
    %
    %   Inputs:
    %       varargin: 
    %           - 1 argument: A vector, which will be converted into 3D array.
    %           - 2 arguments: The first input is a scalar (DC component of phasor array),
    %             and the second is a vector, the positive phasors of the signal. 
    %               The optional argument `z_pospart` is set to `true` by default if not provided.
    %       varg:
    %           - reduce: Boolean flag indicating whether to reduce the array (default: false).
    %           - z_pospart: Boolean flag to enforce positive parts of the phasor (default: false).
    %
    %   Outputs:
    %       obj: A PhasorArray object containing the scalar phasor data.
    %
    %   Notes:
    %       - For a single argument, it must be a vector.
    %       - For two arguments, the first should be a scalar and the second a vector.
    %       - The function defaults to `z_pospart = true` if only two arguments are provided.
    %
    %   Example usage:
    %       obj = ScalarPhasorArray(5, rand(1, 10))  % Creates a phasor array with 5 phasors.
    %
    %   See also: PhasorArray
arguments (Repeating)
    varargin
end
arguments
    varg.reduce=false
    varg.z_pospart=false
end
switch numel(varargin)
    case 1
        assert(isvector(varargin{1}),"for 1 argument, input must be a vector")
        if ~iscolumn(varargin{1})
            varargin{1}=varargin{1}.';
        end
        varargin{1}=permute(varargin{1},[2 3 1]);
    case 2
        assert(isscalar(varargin{1}),"for 2 argument, 1st input must be a scalar")
        assert(isvector(varargin{2}),"for 2 argument, 2nd input must be a vector")
        if ~iscolumn(varargin{2})
            varargin{2}=varargin{2}.';
        end
        varargin{2}=permute(varargin{2},[2 3 1]);
        if ~varg.z_pospart
            warning('two argument provided, switching to z_pos_part=true')
            varg.z_pospart=true;
        end
    otherwise
        error("error check size of input")
end
C=[fieldnames(varg).'; struct2cell(varg).'];
C=C(:).';
obj = PhasorArray(varargin{:},C{:});
end