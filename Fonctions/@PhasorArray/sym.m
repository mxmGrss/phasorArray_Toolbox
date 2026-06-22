function A = sym(nx,ny,h,name,options)
    % SYM Construct a symbolic PhasorArray with structured harmonic components.
    %
    %   A = SYM(NX, NY, H, NAME) creates an NX-by-NY PhasorArray where each
    %   element is a symbolic expression representing harmonic components.
    %
    %   Inputs:
    %     NX   - (integer) Number of rows (default: 1).
    %     NY   - (integer) Number of columns (default: NX).
    %     H    - (integer) Number of harmonics (default: 0).
    %     NAME - (string or cell) Base name for symbolic variables (default: "a").
    %
    %   Outputs:
    %     A - (PhasorArray) Symbolic PhasorArray with structured harmonics.
    %
    %   Notes:
    %     - Each entry in A has a symbolic expression for phasors:
    %       `A_0`, `A_plus_k`, `A_minus_k` for k > 0.
    %     - If NAME is a cell array, it must match the size of NX × NY.
    %     - For scalar NX × NY, NAME is used as the base for all symbols.
    %
    %   See also: sym, ScalarPhasorArray, PhasorArray.
    arguments
        nx=1
        ny=nx
        h = 0
        name ="a"
        options.isreal = false
    end
    %determine if an object is calling or class


    if isa(nx,'PhasorArray')
        A = PhasorArray(sym(nx));
        return
    end
    if (nargin == 1) || (nargin==2 && (ischar(ny) || isstring(ny) || iscell(ny)))
        if ischar(ny) || isstring(ny) || iscell(ny)
            name = ny;
        end
        switch numel(nx)
            case 1
                ny = nx;
                h  = 0;
            case 2
                ny = nx(2);
                nx = nx(1);
                h  = 0;
            case 3
                ny = nx(2);
                h  = nx(3);
                nx = nx(1);
            otherwise
                error('PhasorArray:size:invalidDimInput', 'Dimension argument length must be <=3; got %d.', numel(nx))
        end
    end
    if max(ny,nx)==1
        name={name};
    end
    if iscell(name)
        assert(numel(name)==nx*ny);
        A=PhasorArray.zeros(nx,ny);
        A.Phasor3D = sym(A.Phasor3D);
        for ii =  1:numel(name)
            clear ap am a0 a
            name_i=name{ii};
            ap = sym(name_i+"_plus_",[1 h]);

            if options.isreal
                a0 = sym(name_i+"_0","real");
                a = cat(2,flip(conj(ap)),a0,ap);
            else
                am = sym(name_i+"_minus_",[1 h] );
                a0 = sym(name_i+"_0");
                a = cat(2,flip(am), a0, ap);
            end
            A{ii} = ScalarPhasorArray(a);
        end

    else
        if ~options.isreal
            ap = sym(name+"__%d__%d_plus_%d",[nx ny h]);
            am = sym(name+"__%d__%d_minus_%d",[nx ny h] );
            a0 = sym(name+"__%d__%d_0",[nx ny]);

            a = cat(3,flip(am,3), a0, ap);

            A= PhasorArray(a);
        else
            ap = sym(name+"__%d__%d_plus_%d",[nx ny h]);
            a0 = sym(name+"__%d__%d_0",[nx ny],"real");

            a = cat(3,flip(conj(ap),3), a0, ap);

            A= PhasorArray(a);
        end
    end
end
