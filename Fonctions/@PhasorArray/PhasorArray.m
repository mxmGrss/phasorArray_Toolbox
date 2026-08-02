classdef PhasorArray  < matlab.mixin.indexing.RedefinesParen & matlab.mixin.indexing.RedefinesBrace & matlab.mixin.CustomDisplay
    % PHASORARRAY - Class for handling periodic matrices in the harmonic domain.
    %
    %   The PhasorArray class represents **periodic matrices** via their
    %   **Fourier decomposition**, stored along the third dimension of a 3D array.
    %   The third dimension must be **odd**, with the central slice storing
    %   the **0-th harmonic** component.
    %
    %   This class enables various **mathematical operations, conversions,
    %   time-domain evaluations**, and **harmonic domain analysis**, while
    %   maintaining MATLAB's standard matrix operation compatibility.
    %
    %   This class is used in the PHASORSS class to represent periodic state-space.
    %
    %   --------------------------------------------------------
    %   USAGE EXAMPLES:
    %       A = PhasorArray(ones(3,3,3));     % Create a PhasorArray from a 3D matrix
    %       A.reduce(5);                      % Truncate to 5 harmonics
    %       B = PhasorArray.random(3,3,3);      % Create a random PhasorArray with 3 harmonics (final dimension 3x3x7)
    %       C = A * B;                        % Multiply two PhasorArrays
    %       A.plot();                         % Visualize PhasorArray components
    %       D = C + PhasorArray.random(3,3,4);  % Add two PhasorArrays with different harmonics but compatible first two dimensions
    %   --------------------------------------------------------
    %
    %   DEPENDENCIES:
    %   - Requires MATLAB R2021b or later.
    %   - Uses matlab.mixin.indexing for indexing overloading.
    %
    %   CLASS CONSTRUCTOR:
    %   - obj = PhasorArray(A) creates a PhasorArray from a 3D array A.
    %   - obj = PhasorArray(A, "reduce", true) trims excess zero harmonics.
    %   - obj = PhasorArray(n, m, h) creates an n x m matrix with h harmonics.
    %
    %   --------------------------------------------------------
    %   PROPERTIES:
    %   - Phasor3D: A 3D array representing the harmonics of the PhasorArray.
    %
    %   --------------------------------------------------------
    %   MAIN METHODS:
    %   - **Analysis & Information Retrieval:**
    %       - info (Retrieve PhasorArray properties, dimensions, reality)
    %       - h (Number of harmonics)
    %       - size (Size of PhasorArray)
    %       - isreal (Check if PhasorArray A in time domain A(t) is real)
    %       - dim (Dimension of PhasorArray A in time domain A(t), excluding harmonics)
    %
    %   - **Mathematical Operations:**
    %       - plus, minus, times, rdivide, power, mldivide, mrdivide, mpower
    %       Implement the time domain operation A(t) = B(t) + C(t) for PhasorArrays
    %
    %   - **Harmonic Domain Analysis:**
    %       - Toeplitz Formalism: BT, TB, spBT, spTB
    %       - Fourier Formalism: FvTB, F_bt, F_tb, compatible with Toeplitz
    %       - Eigenvalue Analysis: HmqEig, HmqNEig
    %
    %   - **Time-Domain Evaluation:**
    %       - evalTime, evalp, initial, lsim, sim
    %
    %   - **Phasor Manipulation:**
    %       - reshape, expandBase, extract, pad, reduce, squishBase, trunc
    %
    %   --------------------------------------------------------
    %   COMPLETE METHODS BREAKDOWN:
    %
    %   - **Initialization:**
    %       PhasorArray (Constructor), eye, ones, random,
    %       randomPhasorArrayWithPole, scalar, sym, zeros, empty, ndsdpvar
    %
    %   - **Analysis & Information Retrieval:**
    %       - Element Count and Dimension:
    %           dim, end, h, isempty, length, ndims, numel, numelt, size
    %       - General Analysis:
    %           energy, imagEnergy, realEnergy, pageEnergy, tolReal, info
    %       - Structural Analysis:
    %           iscolumn, iscomplex, ishermitian, isimag, islogical,
    %           ismatrix, isnumeric, isreal, isrow, isscalar, issquare,
    %           issymmetric, isvector, iszero
    %
    %   - **Conversion and Representation:**
    %       - Phasor Representation:
    %           AngleAmpForm, CosSinForm, ImagRealForm, RealImagForm,
    %           SinCosForm
    %       - PhasorArray Special Constructors:
    %           fromTBMatrix, fromF_tb, funcToPhasorArray, time2Phasor,
    %           cqt2ScalarPhasor
    %       - General Conversion and Extraction:
    %           double, sdpval, value, Value, squeeval
    %
    %   - **Harmonic Domain Analysis:**
    %       - Eigenvalue Analysis:
    %           HmqEig, HmqNEig
    %       - Control:
    %           place
    %       - Toeplitz Formalism:
    %           FvTB, F_bt, F_tb, spBT, spBTHankel, spTB, spTBHankel,
    %           spTBmtimes, BT, BTHankel, TB, TBHankel, TBmtimes
    %
    %   - **Mathematical Operations:**
    %       - Arithmetics:
    %           ctroplus, ctrretro, d, det, expm, inv, kron, logm,
    %           mctranspose, mpower, mtranspose, oplus, power, sum,
    %           trace, troplus, ctranspose, transpose, minus, mldivide,
    %           mrdivide, mtimes, plus, rdivide, times, uminus, uplus,
    %           antiD, retro, trretro, ldivide, PhaseShift, pmax
    %       - General:
    %           diag, mconj, mimag, mreal
    %
    %   - **Operator Overloading:**
    %       - Logical:
    %           and, eq, ge, gt, le, lt, ne, not, or
    %
    %   - **Phasor Manipulation:**
    %       - Reshaping:
    %           cat, horzcat, vect, flip, fliplr, flipud, permute,
    %           repeat, repmat, reshape, rot90, vertcat, blkdiag
    %       - Page Operation on Phasor Elements:
    %           pageldivide, pageplus, pagepower, pagerdivide, pagetimes,
    %           pageabs, pageconj, pageimag, pagereal, pagectranspose,
    %           pagetranspose
    %       - Reduction and Expansion:
    %           expandBase, extract, neglect, pad, reduce, squishBase,
    %           trunc
    %       - General:
    %           phas, phasAssign, sub
    %
    %   - **Time-Domain Evaluation:**
    %       evalTime, evalp, initial, lsim, sim
    %
    %   - **Indexing Overloading:**
    %       - Parenthesis:
    %           parenAssign, parenDelete, parenListLength, parenReference
    %       - Brace:
    %           braceAssign, braceListLength, braceReference
    %
    %   - **Utilities:**
    %       info
    %
    %   - **Deprecated:**
    %       confirm_reality
    %
    %   - **Visualization:**
    %       barsurf, plot, plot3D, stem
    %
    %
    %   --------------------------------------------------------
    %   VERSION & AUTHOR:
    %   - Version: 1.0
    %   - Author: Maxime GROSSO
    %   - Contact: maxime.grosso@protonmail.com
    %   --------------------------------------------------------
    %
    %   See also PHASORSS, SPARSEPHASORARRAY
    properties
        Phasor3D %A 3D array reprensenting the harmonics of PhasorArray. Can be a 3D Array of double, of sim or sdpvar
    end
    methods
        function obj = PhasorArray(varargin,nvp)
            % PHASORARRAY - Constructor for the PhasorArray class.
            %
            %   Creates a PhasorArray object representing a periodic matrix using
            %   its Fourier decomposition. The third dimension of the input array
            %   must be **odd**, as it represents harmonics with the **central slice**
            %   storing the **0-th harmonic**.
            %
            %   --------------------------------------------------------
            %   SYNTAX:
            %       pA = PhasorArray(A)
            %           - A is a 3D array (size(A,3) must be odd)
            %           - Supports **double, sym, or sdpvar** data types
            %
            %       pA = PhasorArray(A, "isreal", true)
            %           - A is a 3D array where A(:,:,1) is the zero phasor
            %           - A(:,:,k) contains positive phasors, and negative ones
            %             are inferred via complex conjugation
            %
            %       pA = PhasorArray(A0, Ap, "isreal", true)
            %           - A0 is a 2D array representing the **zero phasor**
            %           - Ap is a 3D array storing **positive harmonics**
            %           - Negative harmonics are inferred via complex conjugation
            %
            %       pA = PhasorArray(A, "reduce", true)
            %           - Automatically removes **zero-valued outer harmonics**
            %
            %       pA = PhasorArray(n, m)
            %           - Creates a **zero-filled** PhasorArray of size `(n, m, 1)`
            %
            %   --------------------------------------------------------
            %   INPUTS:
            %   - varargin (Repeating Arguments):
            %       * A (3D array, size(A,3) must be odd)
            %       * n, m (Scalars specifying dimensions for a zero PhasorArray)
            %       * if using "isreal" flag:
            %           - A0 (2D array): Zero phasor in conjuction with Ap
            %           - Ap (3D array): Positive phasors in conjuction with A0
            %           - A (3D array): Zero and positive phasors at once
            %
            %   - nvp (Name-Value Arguments):
            %       * "reduce" (true/false)
            %           - Removes symmetric zero phasors (default: false)
            %       * "isreal" (true/false)
            %           - Indicates that only zero and positive phasors are provided,
            %             and negative phasors should be inferred (default: false)
            %
            %   --------------------------------------------------------
            %   OUTPUT:
            %   - obj: A `PhasorArray` object storing the Fourier decomposition.
            %
            %   --------------------------------------------------------
            %   EXAMPLES:
            %       % Construct from a random 3D array
            %       A = rand(4,4,5);
            %       pA = PhasorArray(A);
            %
            %       % Construct using zero and positive phasors
            %       A0 = eye(3);
            %       Ap = rand(3,3,2);
            %       pA = PhasorArray(A0, Ap, "isreal", true);
            %
            %       % Create a zero-initialized PhasorArray
            %       pA = PhasorArray(5, 5);
            %
            %   --------------------------------------------------------
            arguments (Repeating)
                varargin
            end
            arguments
                nvp.reduce=false
                nvp.isreal=false
            end
            if nargin == 0
                obj = PhasorArray(1);
                return
            end

            if nvp.isreal % if we indicated to provide a 0-phasor and positive phasors
                if numel(varargin)==1
                    varP=varargin{1};
                    varP0=varP(:,:,1);
                    varPpos=varP(:,:,2:end);
                    obj = PosPart2PhasorArray(varP0,varPpos);
                else
                    if numel(varargin) ~= 2
                        error('PhasorArray:InvalidArgumentCount', 'If 0 and PosPhasor are provided, only 2 arguments must be provided');
                    end
                    obj = PosPart2PhasorArray(varargin{:});
                end
            elseif numel(varargin)>1 % OTHERWISE, if we provided multiple scalar arguments, we're providing a size for a zero phasor array
                if ~(isscalar(varargin{1}) && isscalar(varargin{2}))
                    error('PhasorArray:InvalidArgumentType', 'Arguments for zero matrix must be scalars');
                end
                obj.Phasor3D=zeros(varargin{:});
                nvp.reduce =0;
            elseif isa(varargin{1},'PhasorArray') % if the first argument is a phasor array, all following arguments are ignored
                obj.Phasor3D = pvalue(varargin{1});
            elseif isa(varargin{1},'sparsePhasorArray') % if the first argument is a sparse phasor array, all following arguments are ignored
                obj = varargin{1}.toPhasorArray();
                return
            else
                if mod(size(varargin{1},3),2) ~= 1
                    error('PhasorArray:InvalidDimension', 'Third dimension must be odd (2h+1) to encode h harmonics; got even size %d.', size(varargin{1},3));
                end
                obj.Phasor3D = varargin{1};
            end
            if nvp.reduce
                if isa(obj.Phasor3D,'double')
                    obj.Phasor3D=ReduceArray(obj.Phasor3D);
                end
            end

        end
        function info = info(pA1, tol)
            % INFO Retrieve and evaluate properties of a PhasorArray object
            %   info = INFO(pA1, tol) computes multiple characteristics of the
            %   PhasorArray, including harmonic information, real and imaginary energy
            %   distribution, and shape properties. It returns a structured output
            %   containing these metrics.
            %
            %   INPUTS:
            %       pA1  - The PhasorArray object to analyze.
            %       tol - (Optional) A tolerance used to determine if the object
            %             is considered real. Defaults to machine precision if not
            %             specified.
            %
            %   OUTPUT:
            %       info - A struct containing:
            %          - isReal             : Logical flag indicating if the array
            %                                 is real within the specified tolerance.
            %          - maxImagCoeff       : Maximum magnitude of the imaginary
            %                                 coefficients.
            %          - minImagCoeff       : Minimum magnitude of the imaginary
            %                                 coefficients.
            %          - maxRealCoeff       : Maximum magnitude of the real
            %                                 coefficients.
            %          - minRealCoeff       : Minimum magnitude of the real
            %                                 coefficients.
            %          - iscomplex          : Logical flag indicating if the array
            %                                 has non-zero imaginary content.
            %          - realEnergy_EW      : Element-wise real energy distribution.
            %          - imagEnergy_EW      : Element-wise imaginary energy distribution.
            %          - realEnergy         : Total real energy.
            %          - imagEnergy         : Total imaginary energy.
            %          - Energy_EW          : Element-wise total energy distribution.
            %          - Energy             : Total energy.
            %          - realRelativeEnergy : Fraction of total energy carried by real
            %                                 parts.
            %          - imagRelativeEnergy : Fraction of total energy carried by
            %                                 imaginary parts.
            %          - imagRelativeEnergy_EW : Element-wise relative imaginary energy.
            %          - realRelativeEnergy_EW : Element-wise relative real energy.
            %          - h                  : Number of harmonics in the array.
            %          - size               : Size of the PhasorArray.
            %          - isSymmetric        : Logical flag for symmetry check.
            %          - isHermitian        : Logical flag for Hermitian check.
            %          - isSquare           : Logical flag for square dimension check.
            %          - isScalar           : Logical flag for scalar dimension check.
            %          - isVector           : Logical flag for vector dimension check.
            %
            %   NOTE:
            %       If the array contains small imaginary components but is flagged
            %       as real, proceed with caution if the imaginary energy is non-zero.
            arguments
                pA1
                tol  = []
            end
            if nargin == 1
                warning('PhasorArray:isReal:defaultTol', 'tol not specified, defaulting to eps. Real-valued evaluation may be inaccurate.')
            end
            %display information about the PhasorArray
            fprintf("PhasorArray of size [%d, %d, %d] with %d harmonics\n", size(pA1, 1), size(pA1, 2), size(pA1, 3), pA1.h);
            [r,~] = isreal(pA1,tol);
            info.isReal = r;
            info.tolReal = tolReal(pA1);
            info.tolZero = tolZero(pA1);
            info.maxImagCoeff = max(abs(value(mimag(pA1))),[],'all');
            info.minImagCoeff = min(abs(value(mimag(pA1))),[],'all');
            info.maxRealCoeff = max(abs(value(mreal(pA1))),[],'all');
            info.minRealCoeff = min(abs(value(mreal(pA1))),[],'all');

            info.iscomplex = ~r;

            ro1 = mreal(pA1);
            info.realEnergy_EW = sum(ro1.value.*conj(ro1.value),3);
            info.imagEnergy_EW = sum(value(mimag(pA1)).*conj(value(mimag(pA1))),3);
            info.imagEnergy = sum(info.imagEnergy_EW,'all');
            info.realEnergy = sum(info.realEnergy_EW,'all');
            info.Energy_EW = info.realEnergy_EW + info.imagEnergy_EW;
            info.Energy = sum(info.Energy_EW,'all');
            info.realRelativeEnergy = info.realEnergy/info.Energy;
            info.imagRelativeEnergy = info.imagEnergy/info.Energy;
            info.imagRelativeEnergy_EW = info.imagEnergy_EW./info.Energy_EW;
            info.realRelativeEnergy_EW = info.realEnergy_EW./info.Energy_EW;



            info.h = pA1.h;
            info.size = size(pA1);
            info.isSymmetric = issymmetric(pA1);
            info.isHermitian =ishermitian(pA1);
            info.isSquare = issquare(pA1);
            info.isScalar = isscalar(pA1);
            info.isVector = isvector(pA1);

            if info.imagEnergy>0 && info.isReal
                fprintf("PhasorArray is real within machine precision, but has non zero imaginary energy, proceed with caution, tolReal is %e",tolReal(pA1));
            end

        end
        function [Eew,E] = realEnergy(pA1,elementwise)
            %REALENERGY Compute the total real energy of the PhasorArray.
            %   [Eew,E] = REALENERGY(pA1, elementwise) returns element-wise real energy Eew
            %   and total real energy E.
            %   If elementwise is true, returns element-wise energy. Default is false.
            %
            %   Inputs:
            %       pA1 - The PhasorArray object.
            %       elementwise - (Optional, for single output) Logical flag indicating whether to return element-wise energy.
            %                     If true, returns element-wise energy. Default is false.
            %                     If two outputs are requested, the total energy is returned as second output regardless of this flag.
            %
            %   Outputs:
            %       Eew - Element-wise real energy.
            %       E   - Total real energy.
            %
            %   Example:
            %       pa = PhasorArray(rand(3,3,5));
            %       [Eew, E] = realEnergy(pa);
            %
            %   See also: energy, imagEnergy, DCenergy

            if nargin < 2, elementwise = false; end
            [Eew, E] = energyOf(pA1.mreal, elementwise, nargout);
        end

        function [Eew,E] = imagEnergy(pA1,elementwise)
            %IMAGENERGY Compute the total imaginary energy of the PhasorArray.
            %   [Eew,E] = IMAGENERGY(pA1, elementwise) returns element-wise imaginary energy Eew
            %   and total imaginary energy E.
            %   If elementwise is true, returns element-wise energy. Default is false.
            %
            %   Inputs:
            %       pA1 - The PhasorArray object.
            %       elementwise - (Optional, for single output) Logical flag indicating whether to return element-wise energy.
            %                     If true, returns element-wise energy. Default is false.
            %                     If two outputs are requested, the total energy is returned as second output regardless of this flag.
            %
            %   Outputs:
            %       Eew - Element-wise imaginary energy.
            %       E   - Total imaginary energy.
            %
            %   Example:
            %       pa = PhasorArray(rand(3,3,5));
            %       [Eew, E] = imagEnergy(pa);
            %
            %   See also: energy, realEnergy, DCenergy

            if nargin < 2, elementwise = false; end
            [Eew, E] = energyOf(pA1.mimag, elementwise, nargout);
        end

        function [Eew,E] = energy(pA1,elementwise)
            %ENERGY Compute the total energy of the PhasorArray.
            %   [Eew,E] = ENERGY(pA1, elementwise) returns element-wise energy Eew
            %   and total energy E.
            %   If elementwise is true, returns element-wise energy. Default is false.
            %
            %   The energy is computed as: the \ell_2 norm squared of the PhasorArray over all harmonics.
            %   Through Parseval's theorem, this corresponds to the total energy in the time domain : E = \int |A(t)|^2 dt over one period.
            %
            %   formula E(i,j) = sum |pA1(i,j,k)|^2
            %           E = sum E(i,j)
            %
            %   Inputs:
            %       pA1 - The PhasorArray object.
            %       elementwise - (Optional, for single output) Logical flag indicating whether to return element-wise energy.
            %                     If true, returns element-wise energy. Default is false.
            %                     If two outputs are requested, the total energy is returned as second output regardless of this flag.
            %
            %   Outputs:
            %       Eew - Element-wise total energy.
            %       E   - Total energy.
            %
            %   Example:
            %       pa = PhasorArray(rand(3,3,5));
            %       [Eew, E] = energy(pa);
            %
            %   See also: realEnergy, imagEnergy, DCenergy

            if nargin < 2, elementwise = false; end
            [Eew, E] = energyOf(pA1, elementwise, nargout);
        end
        function [Eew,E] = ACenergy(pA1,elementwise)
            %ACENERGY Compute the energy of the AC (non-DC) components of a PhasorArray.
            %   [Eew, E] = ACENERGY(pA1, elementwise) returns the element-wise energy Eew and total energy E
            %   of the AC components (non-DC harmonics) of the PhasorArray pA1.
            %
            %   Inputs:
            %       pA1 - The PhasorArray object.
            %       elementwise - (Optional, for single output) Logical flag indicating whether to return element-wise energy.
            %                     If true, returns element-wise energy. Default is false.
            %                     If two outputs are requested, the total energy is returned as second output regardless of this flag.
            %
            %   Outputs:
            %       Eew - Element-wise energy of the AC components.
            %       E   - Total energy of the AC components.
            %
            %   Example:
            %       pa = PhasorArray(rand(3,3,5));
            %       [Eew, E] = ACenergy(pa);
            %
            %   See also: energy, realEnergy, imagEnergy, DCenergy
            if nargin < 2, elementwise = false; end
            [Eew, E] = energyOf(pA1.extract(1:pA1.h), elementwise, nargout);
        end
        function [Eew,E] = DCenergy(pA1,elementwise)
            %DCENERGY Compute the energy of the DC (0-th harmonic) component of a PhasorArray.
            %   [Eew, E] = DCENERGY(pA1, elementwise) returns the element-wise energy Eew and total energy E
            %   of the DC component (0-th harmonic) of the PhasorArray pA1.
            %
            %   Inputs:
            %       pA1 - The PhasorArray object.
            %       elementwise - (Optional, for single output) Logical flag indicating whether to return element-wise energy.
            %                     If true, returns element-wise energy. Default is false.
            %                     If two outputs are requested, the total energy is returned as second output regardless of this flag.
            %
            %   Outputs:
            %       Eew - Element-wise energy of the DC component.
            %       E   - Total energy of the DC component.
            %
            %   Example:
            %       pa = PhasorArray(rand(3,3,5));
            %       [Eew, E] = DCenergy(pa);
            %
            %   See also: energy, realEnergy, imagEnergy, ACenergy
            if nargin < 2, elementwise = false; end
            [Eew, E] = energyOf(pA1.extract(0), elementwise, nargout);
        end

        function [Eew, E] = pageEnergy(pA1, nvp)
            % PAGEENERGY Compute and optionally plot the per-harmonic energy of a PhasorArray.
            %
            %   [Eew, E] = pageEnergy(pA1, Name, Value)
            %
            %   Inputs:
            %     pA1           - PhasorArray object
            %
            %   Name-Value Arguments:
            %     'normalized'  - (logical) Normalise outputs by total energy. Default: false.
            %     'excludeDC'   - (logical) Exclude h=0 (DC) from outputs and plot. Default: false.
            %     'cumulative'  - (string)  Cumulative summation mode:
            %                      'none'       — raw per-harmonic values (default)
            %                      'cumulative' — cumsum from h=0 (or h=1 if excludeDC)
            %                      'reverse'    — cumsum from h=H down to 0 (or 1)
            %     'plot'        - (string)  Plot style:
            %                      'none'     — no plot (default)
            %                      'linear'   — line plot, linear scale
            %                      'log'      — line plot, log scale
            %                      'stem'     — stem plot, linear scale
            %                      'stem-log' — stem plot, log scale
            %
            %   Outputs:
            %     Eew  - (rows × cols × N) per-element energy per harmonic, N = h+1 (or h if excludeDC)
            %     E    - (1 × N)           total (Frobenius) energy per harmonic
            %
            %   Example:
            %     [Eew, E] = pageEnergy(A, 'normalized', true, 'excludeDC', true, 'plot', 'stem-log');
            arguments
                pA1
                nvp.normalized  logical = false
                nvp.excludeDC   logical = false
                nvp.cumulative  {mustBeMember(nvp.cumulative,  {'none','cumulative','reverse'})} = 'none'
                nvp.plot        {mustBeMember(nvp.plot, {'none','linear','log','stem','stem-log'})} = 'none'
            end

            % --- 1. Compute per-harmonic energy (one-sided: h=0..H) ---
            for hi = pA1.h : -1 : 0
                [Eew(:,:,hi+1), E(hi+1)] = energy(pA1.extract(hi));
            end

            % --- 2. Exclude DC (h=0) ---
            if nvp.excludeDC
                Eew = Eew(:,:,2:end);   % drop h=0 slice
                E   = E(2:end);
                h_start = 1;
            else
                h_start = 0;
            end

            % --- 3. Normalise ---
            if nvp.normalized
                Etotal   = sum(E(:));
                EewTotal = sum(Eew, 3);
                if Etotal > 0
                    Eew = Eew ./ EewTotal;
                    E   = E   ./ Etotal;
                end
                Eew(isnan(Eew)) = 0;
                E(isnan(E))     = 0;
            end

            % --- 4. Cumulative summation ---
            switch nvp.cumulative
                case 'cumulative'
                    Eew = cumsum(Eew, 3);
                    E   = cumsum(E);
                case 'reverse'
                    Eew = flip(cumsum(flip(Eew, 3), 3), 3);
                    E   = flip(cumsum(flip(E)));
                    % 'none' : nothing to do
            end

            % --- 5. Plot ---
            if ~strcmp(nvp.plot, 'none')
                useLog  = ismember(nvp.plot, {'log',  'stem-log'});
                useStem = ismember(nvp.plot, {'stem', 'stem-log'});

                norm_lbl = ""; if nvp.normalized; norm_lbl = "Normalised "; end
                dc_lbl   = ""; if nvp.excludeDC;  dc_lbl   = " (AC only)"; end

                switch nvp.cumulative
                    case 'cumulative'
                        Eplot = E;
                        ttl   = norm_lbl + "Cum(h) = \Sigma_{k\leq h} E_k" + dc_lbl;
                    case 'reverse'
                        Eplot = E;
                        ttl   = norm_lbl + "ReCum(h) = \Sigma_{k\geq h} E_k" + dc_lbl;
                    otherwise
                        Eplot = E;
                        ttl   = norm_lbl + "E(h) — energy per harmonic" + dc_lbl;
                end

                hx = h_start : h_start + numel(E) - 1;
                if useStem
                    stem(hx, Eplot, 'filled');
                else
                    plot(hx, Eplot);
                end
                title(ttl);
                xlabel("Harmonic order");
                ylabel("Energy");
                if useLog; set(gca, 'YScale', 'log'); end
            end
        end




        function r = plus(pA1,pA2)
            % PLUS Overloads the plus (+) operator for the PhasorArray class.
            %    PLUS(A,B) returns a new PhasorArray representing the element-wise
            %    sum of A(t) and B(t). The second argument can be a scalar or a
            %    scalar PhasorArray, which is broadcast to match the dimensions of A.
            if isscalar(pA2)
                pA2=ones(size(pA1,1),size(pA1,2))*pA2;
            end
            r = PhasorArray(PhasorArrayAdd(pA1,pA2));
            if ~isspecial(r)
                r = r.neglect(0,"exclude0Phasor",false,"reduceMethod","absolute");
            end
        end
        function r = pageplus(pA1,pA2)
            %PAGEPLUS(A,B) : A+B = A(t)+B(t)
            %restricted to phasor of same third dimension
            %shouldn't be used directly, use PLUS instead
            %
            %See also PLUS
            r=pA1.value+pA2.value;
        end
        function r = minus(pA1,pA2)
            %MINUS overloading for PhasorArray
            %MINUS(A,B) : A-B = A(t)-B(t)
            %MINUS can accept scalar or scalar PhasorArray as second argument
            pA2 = -pA2;
            r = pA1 + pA2;
        end
        function r = uminus(pA1)
            %UMINUS overloading for PhasorArray
            %UMINUS(A) : -A = -A(t)
            r = PhasorArray(-pvalue(pA1));
        end
        function r = uplus(pA1)
            %operator overloading for consistency, does nothing
            r = pA1;
        end
        function r = pagetimes(pA1,pA2)
            % Performs element-wise multiplication of phasors from two PhasorArrays
            d=PhasorUnif(pA1,pA2);
            r = PhasorArray(pvalue(d{1}).*pvalue(d{2}));
        end
        function r = times(pA1,pA2)
            %TIMES Element-wise multiplication of two PhasorArray objects.
            %   R = TIMES(O1,O2) returns a PhasorArray containing the element-wise
            %   multiplication of O1 and O2. Both inputs must be PhasorArray objects
            %   (or unify to the same dimension).
            % Assembled with cat rather than written into a preallocated array:
            % zeros() is double, and assigning an sdpvar into it dropped the
            % decision variables silently, a sym one raised outright.
            % PhasorArrayTimes reads each operand's order on its own, so padding
            % them to a common one beforehand only inflated the result -- a
            % constant mask against an order-5 array came back at order 10.
            a = pvalue(PhasorArray(pA1));
            b = pvalue(PhasorArray(pA2));
            [n1, n2, ~] = size(a);
            rows = cell(n1, 1);
            for ii = 1:n1
                cols = cell(1, n2);
                for jj = 1:n2
                    cols{jj} = PhasorArrayTimes(a(ii,jj,:), b(ii,jj,:), [], ...
                        "reduce", false, "output", "Array");
                end
                rows{ii} = cat(2, cols{:});
            end
            r = PhasorArray(cat(1, rows{:}));
            try
                if ~isreal(r) && isreal(pA1) && isreal(pA2) %warning very costly for big phasorArray
                    r = mreal(r);
                end
            catch
            end
        end
        function r = triu(pA1, k)
            %TRIU Keep the upper triangular part of A(t).
            %   TRIU(A) zeroes every entry below the main diagonal, at every
            %   instant. TRIU(A, K) zeroes below the K-th diagonal.
            %
            %   A sparsity pattern is not a symmetry, so it has no name in the
            %   PHASORSYMMETRY vocabulary. It is a product with a constant mask,
            %   which leaves the harmonic order alone and preserves realness. The
            %   payload class survives, so this works on sym and sdpvar too.
            %
            %   Example
            %       U = triu(PhasorArray.random(3, 3, 5));
            %
            %   See also tril, diag, phasorSymmetry.
            arguments
                pA1 PhasorArray
                k (1,1) {mustBeInteger} = 0
            end
            r = pA1 .* PhasorArray(triu(ones(size(pA1,1), size(pA1,2)), k));
        end

        function r = tril(pA1, k)
            %TRIL Keep the lower triangular part of A(t).
            %   TRIL(A) zeroes every entry above the main diagonal, at every
            %   instant. TRIL(A, K) zeroes above the K-th diagonal.
            %
            %   See also triu, diag, phasorSymmetry.
            arguments
                pA1 PhasorArray
                k (1,1) {mustBeInteger} = 0
            end
            r = pA1 .* PhasorArray(tril(ones(size(pA1,1), size(pA1,2)), k));
        end

        function r = mtimes(pA1,pA2)
            %MTIMES Overloads the matrix multiplication operator for PhasorArray.
            %   R = MTIMES(O1,O2) performs the time-domain product of the two
            %   PhasorArray objects O1 and O2 (convolution of the 3D arrays along the third dimension), returning a PhasorArray result.
            try
                r = PhasorArray(PhasorArrayTimes(pA1,pA2));
                try
                    if ~isspecial(r) && ~isreal(r) && isreal(PhasorArray(pA1)) && isreal(PhasorArray(pA2))
                        r = mreal(r);
                    end
                catch
                end
            catch e
                if isa(pA2,"PhasorSS")
                    pA1 = PhasorSS([],[],[],pA1);
                    r = pA1 * pA2;
                else
                    rethrow(e)
                end
            end

        end
        function r = rdivide(pA1,pA2,varargin)
            % RDIVIDE Overloads the right array division operator for the PhasorArray class.
            %   RDIVIDE(A,B) returns a PhasorArray that represents A./B term by term.
            for ii = 1:size(pA2,1)
                for jj = 1:size(pA2,2)
                    d2(ii,jj,:)=PhasorInv(pA2(ii,jj,:),varargin{:});
                end
            end
            d2=PhasorArray(d2);
            d=PhasorUnif(pA1,d2);
            r = PhasorArray(d{1}.*d{2});
        end

        function r = ldivide(pA1,pA2,varargin)
            % LDIVIDE Overloads the left array division operator for the PhasorArray class.
            %   LDIVIDE(A,B) returns a PhasorArray that represents A.\B term by term.
            for ii = 1:size(pA1,1)
                for jj = 1:size(pA1,2)
                    d1(ii,jj,:)=PhasorInv(pA1(ii,jj,:),varargin{:});
                end
            end
            d1=PhasorArray(d1);
            d=PhasorUnif(d1,pA2);
            r = PhasorArray(d{1}.*d{2});
        end
        function r = pagerdivide(pA1,pA2)
            % Performs element-wise right division of phasors from 2 PhasorArrays
            d=PhasorUnif(pA1,pA2);
            r = PhasorArray(d{1}./d{2});
        end
        function r = pageldivide(pA1,pA2)
            % Performs element-wise left division of phasors from 2 PhasorArrays
            d=PhasorUnif(pA1,pA2);
            r = PhasorArray(d{1}.\d{2});
        end
        function r = mrdivide(pA1,pA2,nvp)
            %MRDIVIDE Overloads the right matrix division operator (/) for PhasorArray.
            %   X = B / A solves X(t) * A(t) = B(t) directly in the harmonic domain
            %   using mrHmcDivide with adaptive harmonic selection.
            %
            %   By default, it uses 'autoUpdateh', true.
            %
            %   Optional Name-Value Parameters:
            %       'thresholdResidual' - Target relative residual norm (default: 5e-4).
            %       'autoUpdateh'      - Automatically increase harmonics (default: true).
            %       'maxh'             - Hard cap for auto-update (default: 20*InitialH).
            %       'plotConvergence'  - Show diagnostic plots (default: false).
            %       'verbose'          - Console output level (0, 1, 2).
            %
            %   See also: mrHmcDivide, mlHmcDivide, mldivide
            arguments
                pA1
                pA2
                nvp.h                = []
                nvp.thresholdResidual = 5e-4
                nvp.autoUpdateh      = true
                nvp.maxh             = []
                nvp.stagnationWindow = 15
                nvp.stagnationRatio  = 0.02
                nvp.verbose          (1,1) logical = false
                nvp.updateMethod    {mustBeMember(nvp.updateMethod,{'adaptive','incremental'})} = 'adaptive'
            end
            C = namedargs2cell(nvp);
            [r, residual] = mrHmcDivide(PhasorArray(pA1), PhasorArray(pA2), C{:});
            warnIfNotConverged(residual, 'mrHmcDivide(B, A, ...)');
        end
        function r = mldivide(pA1,pA2,nvp)
            %MLDIVIDE Overloads the left matrix division operator (\) for PhasorArray.
            %   X = A \ B solves A(t) * X(t) = B(t) directly in the harmonic domain
            %   using mlHmcDivide with adaptive harmonic selection.
            %
            %   By default, it uses 'autoUpdateh', true.
            %
            %   Optional Name-Value Parameters:
            %       'thresholdResidual' - Target relative residual norm (default: 5e-4).
            %       'autoUpdateh'      - Automatically increase harmonics (default: true).
            %       'maxh'             - Hard cap for auto-update (default: 20*InitialH).
            %       'plotConvergence'  - Show diagnostic plots (default: false).
            %       'verbose'          - Console output level (0, 1, 2).
            %
            %   See also: mlHmcDivide, mrHmcDivide, mrdivide
            arguments
                pA1
                pA2
                nvp.h                = []
                nvp.thresholdResidual = 5e-4
                nvp.autoUpdateh      = true
                nvp.maxh             = []
                nvp.stagnationWindow = 15
                nvp.stagnationRatio  = 0.02
                nvp.verbose          = false
                nvp.updateMethod    {mustBeMember(nvp.updateMethod,{'adaptive','incremental'})} = 'adaptive'
            end
            C = namedargs2cell(nvp);
            [r, residual] = mlHmcDivide(PhasorArray(pA1), PhasorArray(pA2), C{:});
            % The operator form takes no options, so a caller who hits a
            % non-converged solve has nowhere to turn unless told where.
            warnIfNotConverged(residual, 'mlHmcDivide(A, B, ...)');
        end


    function [r,residual] = mrHmcDivide(B, A, nvp)
        % mrHmcDivide solves X(t)*A(t) = B(t) (equiv. X(t) = B(t) / A(t)) directly in the harmonic domain.
        %   X = mrHmcDivide(B, A) computes the harmonic division using Toeplitz matrices:
        %   F(X)_tb * T(A)_tb = F(B)_tb
        %
        %   Fallback to lrHmcDivide noticing that X*A = B <=> A.'*X.' = B.'
        %
        %   Note: This solves X(t)*A(t) = B(t). It is equivalent to, but much faster than,
        %   lyap(A, 0, -B, "h", h, "T", Inf).
        %
        %   Name-Value arguments (forwarded to mlHmcDivide):
        %       "h"                - Harmonic truncation. Default: max(A.h, B.h)
        %       "thresholdResidual"- Convergence threshold for residual norm. Default: 5e-4
        %       "autoUpdateh"      - Adaptively increase h until convergence. Default: false
        %       "maxh"             - Hard upper bound on h when autoUpdateh=true. Default: h0*20
        %       "stagnationWindow" - Look-back window for stagnation detection. Default: 15
        %       "stagnationRatio"  - Min relative improvement to avoid stagnation flag. Default: 0.02
        %       "updateMethod"     - h update rule: 'adaptive' | 'incremental'
        %       "plotConvergence"  - Plot h/residual convergence diagnostics
        %       "verbose"          - Print convergence info. Default: 1

        arguments
            B PhasorArray
            A PhasorArray
            nvp.h                = []
            nvp.thresholdResidual = 5e-4
            nvp.autoUpdateh      = true   % as mlHmcDivide: the mirror operation
                                          % had the opposite default for no reason
            nvp.maxh             = []    % hard upper bound on h (default: h0 * 20)
            nvp.stagnationWindow = 15     % look-back window for stagnation detection
            nvp.stagnationRatio  = 0.02  % relative improvement threshold (< 5% = stagnation)
            nvp.verbose          (1,1) logical = false
            nvp.updateMethod    {mustBeMember(nvp.updateMethod,{'adaptive','incremental'})} = 'adaptive'
        end

        % Same refusal mldivide already gets from sparray2TBlocks, said here so the
        % right-hand division does not answer with a MATLAB internal message.
        solvedPayload(pvalue(A), "mrHmcDivide");
        solvedPayload(pvalue(B), "mrHmcDivide");
        C = namedargs2cell(nvp);

        [r,residual] = mlHmcDivide(A.', B.', C{:});
        r = r.';
    end

    function F = spF_tb(obj, h)
        % spF_tb Returns the sparse Fourier block column matrix F_tb(X, h).
        if nargin < 2 || isempty(h)
            h = obj.h;
        end
        F = sparse(obj.F_tb(h));
    end

    function F = spF_bt(obj, h)
        % spF_bt Returns the sparse Fourier block row matrix F_bt(X, h).
        if nargin < 2 || isempty(h)
            h = obj.h;
        end
        F = sparse(obj.F_bt(h));
    end

    function r = pagepower(pA1,m)
        % pagepower Element-wise exponentiation of phasors.
        %
        %   r = pagepower(pA1, m) computes `pA1 .^ m` phasor-wise, applying exponentiation
        %   to each phasor coefficient individually.
        %
        %   Inputs:
        %   - pA1: A PhasorArray.
        %   - m : A scalar or a PhasorArray of matching size.
        %
        %   Behavior:
        %   - If `m` is a scalar, applies `.^m` to each individual phasor of `pA1`.
        %   - If `m` is a PhasorArray, exponentiation is applied to matching phasor coefficients.
        %   - Each frequency component of `pA1` is exponentiated separately.
        %
        %   Notes:
        %   - This is **not** equivalent to exponentiating the time-domain matrix \( A(t) \).
        %   - Rarely useful in control or signal processing applications.
        %
        %   See also: power, mpower
        if isa(m,'PhasorArray')
            d=PhasorUnif(pA1,m);
            m=d{2};
            r = PhasorArray((pvalue(d{1})).^(pvalue(m)));
        else
            r = PhasorArray(pvalue(pA1).^m);
        end
    end
    function r = power(pA1,m)
        % POWER Element-wise power of a periodic matrix A(t).
        %
        %   r = power(pA1, m) computes the element-wise power of each matrix entry of a periodic matrix \( A(t) \),
        %   represented as a PhasorArray. This is equivalent to computing \( A(t)^{m} \) for each individual element of the matrix
        %   at each time instance in the time domain.
        %
        %   Inputs:
        %   - pA1: A PhasorArray representing the periodic matrix \( A(t) \).
        %   - m: The exponent (scalar value) to apply element-wise.
        %
        %   Output:
        %   - r: A PhasorArray representing the matrix \( A(t)^{m} \) computed element-wise.
        %
        %   Behavior:
        %   - The function applies the power operation **element-wise** for each entry in the time-domain representation of \( A(t) \).
        %   - This corresponds to raising each individual element \( A_{ij}(t) \) of the matrix to the power \( m \).
        %
        %   See also: PhasorPow, mpower, expm
        r=pA1;
        % element wise power
        for ii=1:size(pA1,1)
            for jj=1:size(pA1,2)
                r{ii,jj}=pA1{ii,jj}^m;
            end
        end
    end
    function r = mpower(pA1,m)
        % MPOWER Matrix power of a periodic matrix A(t), computed for integer exponents.
        %
        %   r = mpower(pA1, m) computes the matrix power \( A(t)^m \) for a periodic matrix \( A(t) \),
        %   represented as a PhasorArray. This operation is only valid for **integer exponents** and is done in the time domain.
        %   It computes the repeated matrix multiplication of \( A(t) \) with itself \( m \) times.
        %
        %   Inputs:
        %   - pA1: A PhasorArray representing the periodic matrix \( A(t) \).
        %   - m: A positive or negative integer exponent representing the power.
        %
        %   Output:
        %   - r: A PhasorArray representing \( A(t)^m \), the matrix raised to the \( m \)-th power.
        %
        %   Behavior:
        %   - For positive integer \( m \), it multiplies \( A(t) \) with itself \( m \) times.
        %   - For negative integer \( m \), it computes the inverse of \( A(t) \) raised to the power \( |m| \) and then multiplies.
        %   - For non integer \(m \), it calls PhasorPow, which perform time domain matrix power and then IFFT.
        %
        %   Notes:
        %   - This function applies matrix exponentiation, unlike `power`, which operates element-wise.
        %
        %   See also: PhasorPow, power, expm
        if mod(m,1)==0 && m>1
            prec=pvalue(pA1);
            for ii=2:m
                prec=PhasorArrayTimes(prec,pA1);
            end
        elseif mod(m,1)==0 && m<0
            d1=PhasorInv(pA1);
            prec=d1;
            for ii=-2:-1:m
                prec=PhasorArrayTimes(prec,d1);
            end
        else
            prec=PhasorPow(pA1,m);
        end
        r=PhasorArray(prec);
    end

    function r = oplus(pA1,pA2)
        %OPLUS Kronecker sum of two PhasorArray objects.
        %   R = OPLUS(O1,O2) returns a PhasorArray representing the Kronecker
        %   sum of O1 and O2.
        %
        %   The kronecker sum A(t) oplus B(t) is usually defined as
        %          A(t) otimes I + I otimes B(t)
        %   where I is the identity matrix of the same size as A(t) and B(t)
        %   and otimes is the usual kronecker product.
        %
        %   see also KRON, TROPLUS, CTROPLUS
        r=PhasorArray(PhasorArrayOplus(pA1,pA2));
    end
    function r = troplus(pA1,pA2)
        %TROPLUS Kronecker sum of the transpose of two PhasorArray objects.
        %   R = TROPLUS(O1,O2) returns a PhasorArray representing the Kronecker
        %   sum of the transpose of O1 and O2.
        %
        %   see also KRON, OPLUS, CTROPLUS
        r=PhasorArray(PhasorArrayOplus(pagetranspose(pvalue(pA1)),pA2));
    end
    function r = ctroplus(pA1,pA2)
        %CTROPLUS Kronecker sum of the conjugate transpose of two PhasorArray objects.
        %   R = CTROPLUS(O1,O2) returns a PhasorArray representing the Kronecker
        %   sum of the conjugate transpose of O1 and O2.
        %
        %   see also KRON, OPLUS, TROPLUS
        r=PhasorArray(PhasorArrayOplus(pagectranspose(pvalue(pA1)),pA2));
    end
    function r = kron(pA1,pA2)
        %KRON Kronecker product of two PhasorArray objects.
        %   R = KRON(O1,O2) returns a PhasorArray representing the Kronecker
        %   product of O1 and O2.
        %
        %   see also OPLUS, TROPLUS, CTROPLUS
        r=PhasorArray(PhasorArrayKron(pA1,pA2));
    end
    function r = retro(pA1)
        %RETRO Reverse the time axis of the PhasorArray.
        %   B = RETRO(A) returns a PhasorArray B such that B(t) = A(-t).
        %   This function flips the PhasorArray along the third dimension,
        %   effectively reversing the time axis.
        r=PhasorArray(flip((pvalue(pA1)),3));
    end
    function r = trretro(pA1)
        %TRRETRO Reverse the time axis and transpose the PhasorArray.
        %   B = TRRETRO(A) returns a PhasorArray B such that B(t) = A(-t).'
        %   This function flips the PhasorArray along the third dimension,
        %   effectively reversing the time axis and transposing each page.
        r=PhasorArray(flip(pagetranspose(pvalue(pA1)),3));
    end
    function r = ctrretro(pA1)
        %CTRRETRO Reverse the time axis and conjugate transpose the PhasorArray.
        %   B = CTRRETRO(A) returns a PhasorArray B such that B(t) = A(-t)'.
        %   Conjugating the pages already reverses the harmonic index, so no flip:
        %   with the flip this returned A(t)', duplicating MCTRANSPOSE.
        r=PhasorArray(pagectranspose(pvalue(pA1)));
    end

    function out = extract(pA1,index,symmetric)
        %EXTRACT Extract specified phasors from a PhasorArray.
        %   OUT = EXTRACT(O1, INDEX) returns a PhasorArray containing only the
        %   phasors specified by INDEX. All other phasors are set to zero.
        %
        %   OUT = EXTRACT(O1, INDEX, SYMMETRIC) allows for symmetric extraction.
        %   If INDEX contains only positive integers and SYMMETRIC is true, the
        %   corresponding negative phasors are also included. If INDEX contains
        %   negative integers, SYMMETRIC is forced to false.
        %
        %   INPUTS:
        %       O1        - The PhasorArray object to extract from.
        %       INDEX     - A logical array matching the third dimension of the
        %                   PhasorArray or a vector of positive and negative integers.
        %       SYMMETRIC - (Optional) A logical flag indicating whether to include
        %                   symmetric phasors. Defaults to true.
        %
        %   OUTPUT:
        %       OUT       - A PhasorArray containing only the specified phasors.
        %
        %   EXAMPLES:
        %       A = PhasorArray.random(3, 3, 5);
        %       B = A.extract([1, 3, 5]);
        %       C = A.extract([1, 3], false);
        %       D = A.extract([1, 3], true);
        %
        %   See also trunc, reduce, neglect
        arguments
            pA1
            index
            symmetric=true
        end
        if islogical(index)
            if numel(index)~=size(pA1,3)
                error('PhasorArray:neglect:indexSizeMismatch', 'Logical index size (%d) must match size(obj,3)=%d.', numel(index), size(pA1,3))
            end
            index=find(index);
            %find complement of index
            indexc=setdiff(1:size(pA1,3),index);

            out = pA1;
            out(:,:,indexc)=0;
        else
            if all(index>=0) && symmetric
                index=[index -index];
                maxIndex = max(max(index),pA1.h);
                indexc=setdiff(-maxIndex:maxIndex,index);
                out = pA1;
                out{:,:,indexc}=0;
            else
                out = pA1;
                maxIndex = max(max(abs(index)),pA1.h);
                indexc=setdiff(-maxIndex:maxIndex,index);
                out{:,:,indexc}=0;
            end
        end
    end

    function r = trunc(pA1,m)
        %TRUNC Truncate the PhasorArray to a specified number of phasors.
        %   R = TRUNC(O1, M) returns a new PhasorArray that is truncated to
        %   contain only M phasors. If M is not specified, the function will
        %   use a default value.
        %
        %   INPUTS:
        %       pA1 - The original PhasorArray object to be truncated.
        %       m  - (Optional) The number of phasors to retain in the truncated
        %            PhasorArray. If not provided, a default value is used.
        %
        %   OUTPUT:
        %       r  - The truncated PhasorArray object containing only M phasors.
        %
        %   EXAMPLE:
        %       A = PhasorArray.random(3, 3, 5);
        %       B = A.trunc(3);
        %       % B now contains only 3 phasors from the original PhasorArray A.
        %
        %   See also REDUCE, NEGLECT

        arguments
            pA1
            m=[]
        end
        r=PhasorArray(ReduceArray(pA1,m));
    end
    function r = reduce(pA1,htrunc,nvp)
        %REDUCE Truncate or filter the PhasorArray based on harmonic order or magnitude thresholds.
        %
        %   This function reduces a PhasorArray by either:
        %   - Truncating harmonics beyond a given order (`htrunc`).
        %   - Filtering harmonics based on their magnitude (`reduceMethod` and `reduceThreshold`).
        %   - Optionally controlling how the 0th harmonic is handled (`exclude0Phasor`).
        %   - Applying a hard or soft thresholding strategy (`hardThresholdPhasors`).
        %
        %   Syntax:
        %   r = REDUCE(pA1)
        %       Reduces `pA1` using the default relative thresholding method,
        %       preserving all harmonics unless their magnitudes are below 1e-15.
        %
        %   r = REDUCE(pA1, htrunc)
        %       Truncates `pA1` to harmonics of order ≤ `htrunc`.
        %
        %   r = REDUCE(pA1, htrunc, 'reduceMethod', method, 'reduceThreshold', threshold,
        %                'exclude0Phasor', exclude, 'hardThresholdPhasors', hardThreshold)
        %       Applies a combination of truncation and magnitude filtering based on additional parameters.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The PhasorArray object to be reduced.
        %   - htrunc (integer, optional) : The maximum harmonic order to retain. Default is empty ([]),
        %       meaning harmonics are only removed based on filtering criteria.
        %
        %   Name-Value Pair Arguments:
        %   - 'reduceMethod' (char, optional) : Defines the reduction approach.
        %       - 'absolute' : Removes phasors whose magnitude is below the **absolute** threshold.
        %       - 'relative' (default) : Removes phasors whose magnitude is below a **relative** threshold.
        %
        %   - 'reduceThreshold' (double, optional) : The threshold value for filtering phasors.
        %       - If 'absolute', removes phasors with magnitude < `reduceThreshold`.
        %       - If 'relative', removes phasors with magnitude < `reduceThreshold * max(magnitude)`.
        %       - Default: 1e-15.
        %
        %   - 'exclude0Phasor' (logical, optional) :
        %       - true (default) : The 0th harmonic (DC component) is always considered in max magnitude evaluation.
        %       - false : The 0th harmonic is ignored when computing the relative maximum magnitude.
        %                This affects thresholding in 'relative' mode.
        %
        %   - 'hardThresholdPhasors' (logical, optional) :
        %       - false (default) : `pA1` is truncated to order `m`, meaning **only harmonics beyond order m**
        %         are removed when all phasors of order > `m` fall below the threshold.
        %         Lower-order phasors (< `m`) are **never** removed, even if below threshold.
        %       - true : Any phasor below the threshold is explicitly **set to zero**. Then, trailing zero phasors
        %         (i.e., those bearing no information) are removed to optimize storage.
        %
        %   Output:
        %   - r (PhasorArray) : The reduced PhasorArray after applying truncation or filtering.
        %
        %   Example:
        %   % Truncate PhasorArray to the 5th harmonic
        %   r = reduce(pA1, 5);
        %
        %   % Apply absolute thresholding without truncating harmonics
        %   r = reduce(pA1, [], 'reduceMethod', 'absolute', 'reduceThreshold', 1e-10, 'exclude0Phasor', false);
        %
        %   % Apply soft thresholding with truncation logic
        %   r = reduce(pA1, [], 'hardThresholdPhasors', false, 'reduceMethod', 'relative', 'reduceThreshold', 1e-12);
        %
        %   % Apply hard thresholding where all small phasors are explicitly set to zero
        %   r = reduce(pA1, [], 'hardThresholdPhasors', true, 'reduceMethod', 'relative', 'reduceThreshold', 1e-12);
        %
        %   See also: ReduceArray
        arguments
            pA1
            htrunc=[]
            nvp.reduceMethod char {mustBeMember(nvp.reduceMethod,{'absolute','relative'})} = 'relative'
            nvp.reduceThreshold {mustBeNumeric,mustBeReal} = 1e-15
            nvp.exclude0Phasor (1,1) logical = true
            nvp.hardThresholdPhasors=false
        end
        % Reducing compares magnitudes, which a symbolic or optimisation
        % variable does not have. Left alone the call surfaces a YALMIP
        % strict-inequality error or a MATLAB dimension error, neither of
        % which names the actual problem.
        if ~isnumeric(pvalue(pA1))
            error('PhasorArray:reduce:symbolicPayload', ...
                ['reduce compares phasor magnitudes and needs a numeric array; this one holds %s. ' ...
                 'Build the variable at the order you want, or solve first and reduce the ' ...
                 'numeric result (sdpval for YALMIP, double for symbolic).'], class(pvalue(pA1)));
        end
        r=PhasorArray(ReduceArray(pA1,htrunc,"reduceMethod", nvp.reduceMethod,"reduceThreshold", nvp.reduceThreshold,"exclude0Phasor", nvp.exclude0Phasor,"hardThresholdPhasors", nvp.hardThresholdPhasors));
    end
    function r = neglect(pA1, reduceThreshold, nvp)
        %NEGLECT Set to zero all phasors below a given threshold in a PhasorArray.
        %
        %   This function filters out low-magnitude phasors in a PhasorArray `pA1`,
        %   setting them to zero based on a specified threshold and reduction method.
        %   for matricial phasor array, it acts component-wise, meaning each entry (i,j) is "neglected" independently of the others.
        %
        %   Syntax:
        %   r = NEGLECT(pA1, reduceThreshold)
        %       Sets to zero all phasors in `pA1` with magnitude below `reduceThreshold`,
        %       using the default 'relative' thresholding method.
        %
        %   r = NEGLECT(pA1, reduceThreshold, 'reduceMethod', method, 'exclude0Phasor', exclude, 'h', h)
        %       Applies phasor filtering with additional options.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The PhasorArray object whose small phasors should be neglected.
        %   - reduceThreshold (double) : The magnitude threshold below which phasors are set to zero.
        %
        %   Name-Value Pair Arguments:
        %   - 'reduceMethod' (char, optional) : Defines how the threshold is applied.
        %       - 'absolute' : A phasor is set to zero if its magnitude is below `reduceThreshold`.
        %       - 'relative' (default) : A phasor is set to zero if its magnitude is below
        %         `reduceThreshold * max(magnitude)` (where the reference maximum can exclude the 0th phasor).
        %
        %   - 'exclude0Phasor' (logical, optional) :
        %       - false (default) : The 0th phasor (DC component) **is included** when computing the reference magnitude.
        %       - true : The 0th phasor is **ignored** when computing the reference magnitude for relative thresholding.
        %
        %   - 'h' (integer, optional) :
        %       - If provided, reduces the PhasorArray to at most `h` harmonics after applying the threshold.
        %       - Similar to `reduce(h)`, but performed after zeroing small phasors.
        %
        %   Output:
        %   - r (PhasorArray) : The filtered PhasorArray, with phasors below the threshold set to zero.
        %
        %   Example:
        %   % Set to zero all phasors with magnitude below 1e-15 using default relative method
        %   r = neglect(pA1, 1e-15);
        %
        %   % Apply absolute thresholding (removes all phasors with magnitude < 1e-15)
        %   r = neglect(pA1, 1e-15, 'reduceMethod', 'absolute');
        %
        %   % Use relative thresholding, ignoring the 0th phasor when computing the max reference
        %   r = neglect(pA1, 1e-15, 'reduceMethod', 'relative', 'exclude0Phasor', true);
        %
        %   % Apply thresholding and then truncate to at most 5 harmonics
        %   r = neglect(pA1, 1e-12, 'h', 5);
        %
        %   See also: REDUCE, TRUNC, ReduceArray
        arguments
            pA1
            reduceThreshold {mustBeNumeric,mustBeReal} = 1e-15
            nvp.reduceMethod {mustBeMember(nvp.reduceMethod,{'absolute','relative'})}  = 'relative'
            nvp.exclude0Phasor (1,1) logical = false
            nvp.h=[]
        end
        val=pA1.value;
        h=pA1.h;
        switch  nvp.reduceMethod
            case 'absolute'
                val_rel=val;
            otherwise
                if nvp.exclude0Phasor
                    ref=max(abs(val(:,:,[1:h , (h+2):(2*h+1)])),[],3); %maximum harmonic on each coeef, excepting the phasor 0.
                else
                    ref=max(abs(val),[],3); %maximum harmonic on each coeef, excepting the phasor 0.
                end
                val_rel=abs(val./ref);
        end
        val(abs(val_rel)<reduceThreshold)=0;
        r = PhasorArray(val);
        r = r.reduce(nvp.h,"reduceMethod",nvp.reduceMethod,"reduceThreshold",reduceThreshold,"exclude0Phasor",nvp.exclude0Phasor);
    end


    function r = flip(pA1,dim)
        %FLIP the PhasorArray along the dim dimension
        %   - dim : dimension to flip
        %   - r = flip(pA1,dim) flip the PhasorArray along the dim dimension
        %
        %   See also FLIPLR, FLIPUD
        %
        %   Inputs:
        %       pA1  - PhasorArray object to be flipped
        %       dim - Dimension along which to flip the PhasorArray
        %
        %   Outputs:
        %       r   - Resulting PhasorArray after flipping
        %
        %   Example:
        %       pa = PhasorArray(rand(3,3,3));
        %       pa_flipped = flip(pa, 2);
        %
        %   Note:
        %       If dim is 3, a warning is issued indicating that the flip
        %       along the third dimension produces M(-t).
        if dim==3
            warning('PhasorArray:flip:reversesTime', 'Flipping along dim 3 produces M(-t), i.e., a time-reversal.')
        end
        r=PhasorArray(flip(pA1.value,dim));
    end
    function r = fliplr(pA1)
        %FLIPLR Flip the PhasorArray left to right
        %   r = FLIPLR(pA1) flips the PhasorArray pA1 in the left-right direction.
        %   The function returns a new PhasorArray r with the elements flipped.
        %
        %   Example:
        %       pa = PhasorArray(rand(3,3,3));
        %       pa_flipped = fliplr(pa);
        %
        %   See also FLIP, FLIPUD

        r=PhasorArray(fliplr(pA1.value));
    end
    function r = flipud(pA1)
        %FLIPUD Flip the PhasorArray up to down
        %   r = FLIPUD(pA1) flips the PhasorArray pA1 in the up-down direction.
        %   The function returns a new PhasorArray r with the elements flipped.
        %
        %   Example:
        %       pa = PhasorArray(rand(3,3,3));
        %       pa_flipped = flipud(pa);
        %
        %   See also FLIP, FLIPLR

        r=PhasorArray(flipud(pA1.value));
    end

    function r = rot90(pA1,varargin)
        %ROT90 Rotate the PhasorArray by 90 degrees
        %   r = ROT90(pA1) rotates the PhasorArray object pA1 by 90 degrees
        %   counterclockwise and returns the result in r.
        %
        %   r = ROT90(pA1, k) rotates the PhasorArray object pA1 by 90 degrees
        %   counterclockwise k times.
        %
        %   Input:
        %       pA1 - PhasorArray object to be rotated
        %       k  - (Optional) Number of times to rotate by 90 degrees
        %
        %   Output:
        %       r  - PhasorArray object after rotation
        %
        %   see also FLIP, FLIPLR, FLIPUD
        r=PhasorArray(rot90(pA1.value,varargin{:}));
    end


    function r = horzcat(pA1,varargin)
        %HORZCAT Concatenate PhasorArray objects horizontally
        %
        %   R = HORZCAT(O1, VARARGIN) concatenates the PhasorArray object O1
        %   with additional PhasorArray objects provided in VARARGIN,
        %   horizontally, if their vertical dimensions are compatible.
        %
        %   Input arguments:
        %   O1 - The first PhasorArray object.
        %   VARARGIN - Additional PhasorArray objects to concatenate with O1.
        %
        %   Output arguments:
        %   R - The resulting PhasorArray object after horizontal concatenation.
        %
        %   If only one input argument is provided, the function returns O1.
        %   If the third dimension sizes of O1 and the current object in VARARGIN
        %   are not equal, the function uses PhasorUnif to unify them before
        %   concatenation.
        %
        %   Example:
        %       r = horzcat(pA1, pA2, pA3);
        %
        %   See also: VERTCAT, PHASORUNIF, PVALUE
        if nargin==1
            r=pA1;
            return
        end
        for ii=1:numel(varargin)
            pA2=varargin{ii};
            if size(pA1,3) ~= size(pA2,3)
                d=PhasorUnif(pA1,pA2);
                r = PhasorArray([pvalue(d{1}) pvalue(d{2})]);
            else
                u1=pvalue(pA1);
                u2=pvalue(pA2);
                r = PhasorArray([u1 , u2]);
            end
            pA1=r;
        end
    end
    function r = vertcat(varargin)
        %VERTCAT Concatenate PhasorArray objects vertically
        %   R = VERTCAT(O1, O2, ...) concatenates the PhasorArray object O1
        %   with additional PhasorArray objects provided in VARARGIN,
        %   vertically, if their horizontal dimensions are compatible.
        %
        %   Input arguments:
        %   O1 - The first PhasorArray object.
        %   VARARGIN - Additional PhasorArray objects to concatenate with O1.
        %
        %   Output arguments:
        %   R - The resulting PhasorArray object after vertical concatenation.
        %
        %   If only one input argument is provided, the function returns O1.
        %   If the third dimension sizes of O1 and the current object in VARARGIN
        %   are not equal, the function uses PhasorUnif to unify them before
        %   concatenation.
        %
        %   Example:
        %       r = vertcat(pA1, pA2, pA3);
        %
        %   See also: HORZCAT, PHASORUNIF, PVALUE

        if nargin==1
            r=varargin{1};
            return
        end

        if nargin>=2
            r=cat(1,varargin{:});
            return
        else
            pA1=varargin{1};
            pA2=varargin{2};
        end
        if size(pA1,3)~=size(pA2,3)
            d=PhasorUnif(pA1,pA2);
            r = PhasorArray([pvalue(d{1}) ; pvalue(d{2})]);
        else
            u1=pvalue(pA1);
            u2=pvalue(pA2);
            r = PhasorArray([u1 ; u2]);
        end

    end
    function r = blkdiag(pA1)
        %BLKDIAG of phasorarray
        %  - r = blkdiag(pA1,pA2,...) : block diagonal of phasorArray
        %
        %   See also HORZCAT, VERTCAT
        %
        %   This function creates a block diagonal matrix from the input
        %   phasor arrays. It takes multiple phasor arrays as input and
        %   constructs a block diagonal matrix where each block corresponds
        %   to one of the input phasor arrays.
        %
        %   Input:
        %       pA1 - Repeating argument list of phasor arrays
        %
        %   Output:
        %       r - Resulting PhasorArray object containing the block
        %           diagonal matrix
        %
        %   Example:
        %       pa1 = PhasorArray(rand(3,3,2));
        %       pa2 = PhasorArray(rand(2,2,2));
        %       result = blkdiag(pa1, pa2);
        %
        %   This example creates a block diagonal matrix from two phasor
        %   arrays pa1 and pa2.
        arguments (Repeating)
            pA1
        end
        b1=PhasorUnif(pA1{:});
        h=size(b1{1},3);
        out=cell(h,1);
        for hi=1:h
            argi=cell(numel(b1),1);
            for b1i=1:numel(b1)
                bb1i=pvalue(b1{b1i});
                argi{b1i}=bb1i(:,:,hi);
            end
            out{hi}=blkdiag(argi{:});
        end
        out=cat(3,out{:});
        r=PhasorArray(out);
    end
    function r = repmat(pA1,M,N)
        %REPMAT Replicate and tile a PhasorArray.
        %   r = REPMAT(pA1,M,N) replicates and tiles the PhasorArray pA1 M times along the first dimension and N times along the second dimension.
        %   r = REPMAT(pA1,M) replicates and tiles the PhasorArray pA1 M times along the first and second dimensions.
        %   r = REPMAT(pA1,M,N,L) replicates and tiles the PhasorArray pA1 M times along the first dimension, N times along the second dimension, and L times along the third dimension.
        %
        %   Inputs:
        %       pA1 - The PhasorArray object to be replicated.
        %       M  - Number of times to replicate along the first dimension.
        %       N  - (Optional) Number of times to replicate along the second dimension.
        %
        %   Outputs:
        %       r  - The resulting PhasorArray after replication.
        %
        %   Example:
        %       pa = PhasorArray(rand(3,3,3));
        %       pa_replicated = repmat(pa, 2, 3);
        %
        %   See also REPMAT, RESHAPE

        arguments
            pA1
            M
        end
        arguments (Repeating)
            N
        end
        switch nargin
            case 1
                r = pA1;
            case 2
                switch numel(M)
                    case 1
                        N=M;
                        L=1;
                        M=[M N L];
                    case 2
                        N=M(2);
                        M=M(1);
                        L=1;
                        M=[M N L];
                    case 3
                        warning('PhasorArray:repmat:extraDimIgnored', 'repmat: extra dimension(s) beyond dim 1-2 are not meaningful for PhasorArray. Proceeding with caution.')
                    otherwise
                        warning('PhasorArray:repmat:extraDimIgnored', 'repmat: extra dimension(s) beyond dim 1-2 are not meaningful for PhasorArray. Proceeding with caution.')
                end
                r = PhasorArray(repmat(pA1.value,M));
                return
            case 3
                assert(numel(M)==1)
                N=N{1};
                assert(numel(N)==1)
                L=1;
                r = PhasorArray(repmat(pA1.value,M,N,L));
                return
            otherwise
                warning('PhasorArray:repmat:extraDimIgnored', 'repmat: extra dimension(s) beyond dim 1-2 are not meaningful for PhasorArray. Proceeding with caution.')
                r = PhasorArray(repmat(pA1.value,M,N{:}));
                return
        end

    end
    function r = reshape(pA1,varargin)
        % reshape Reshape a PhasorArray while preserving the phasor dimension.
        %
        %   r = reshape(pA1, M, N)     - Reshapes `pA1` into an MxN PhasorArray,
        %                               preserving the third (phasor) dimension.
        %   r = reshape(pA1, M, N, L)  - Reshapes `pA1` into an MxNxL PhasorArray.
        %
        %   Notes:
        %   - If only `M, N` are provided, `L` is automatically set to `size(pA1,3)`.
        %   - The third dimension (phasor order) should remain unchanged in most cases.
        %
        %   See also: repmat, reshape
        if numel(varargin)==2
            for dimi = 3:ndims(pA1)
                varargin{dimi} = size(pA1,dimi);
            end
        end
        r=PhasorArray(reshape(pvalue(pA1),varargin{:}));
    end
    function r = permute(pA1,varargin)
        % permute Reorder dimensions of a PhasorArray.
        %
        %   r = permute(pA1, order)   - Rearranges the dimensions of `pA1` according to `order`.
        %   r = permute(pA1)          - Swaps the first and second dimensions (acts as transpose).
        %
        %   Notes:
        %   - By default, swaps dimensions `[1 2]`, equivalent to a transpose.
        %   - The third dimension (phasor order) is **preserved** unless explicitly modified.
        %   - Usually applied to **matrix-like** PhasorArrays, not 3D arrays.
        %
        %   See also: transpose, permute
        if nargin==1
            varargin{1}=[2 1];
        end
        if numel(varargin{1})==2
            varargin{1}=[varargin{1} 3];
        end
        r=PhasorArray(permute(pvalue(pA1),varargin{:}));
    end
    function r = sub(pA1,n1,n2)
        % sub Extract specific elements from a PhasorArray based on indices.
        %
        %   r = sub(pA1, n1, n2) extracts the phasor array component
        %   `Phasor3D(n1, n2, :)`, equivalent to `A{n1, n2}`.
        %
        %   Inputs:
        %       pA1  - The PhasorArray object to extract from.
        %       n1  - Row indices or logical mask (can be `:` for all rows).
        %       n2  - (Optional) Column indices (default: 1 if pA1 is a column vector).
        %
        %   Output:
        %       r   - Extracted PhasorArray elements with `reduce=false` to
        %             preserve the structure.
        %
        %   Behavior:
        %   - If `n1` is a logical array, it is converted to indices.
        %   - If `pA1` is a column vector (`size(pA1,2)==1`), `n2` defaults to 1.
        %   - If `n1` and `n2` are omitted, the function reshapes `r` to match
        %     `size(n1,1) × size(n1,2)`.
        %
        %   See also: PhasorArray, reduce, logical indexing
        arguments
            pA1
            n1
            n2=[]
        end
        if nargin<3
            % case where we apparently have only one coordinate entry
            % case 1 it's an integer matrix
            % case 2 it's a logical matrix, and we need to verify
            % that it has the right size

            if ~strcmp(n1,':')% matrix input
                if islogical(n1)
                    if numelt(pA1)==numel(n1)
                        n1=find(logical(n1));
                    else
                        error('PhasorArray:subsasgn:indexSizeMismatch', 'Logical index size %dx%d must match first two dims of PhasorArray %dx%d.', size(n1,1), size(n1,2), size(pA1,1), size(pA1,2))
                    end
                end
            end

            if size(pA1,2)==1 % its already a column PhasorArray, normal to have only one index
                n2=1;
            else
                pA1=vect(pA1); % we make it a column, it will be easier to manipulate
                n2=1;
            end
        end
        indices = repmat({':'}, 1, ndims(pA1));
        indices{1} = n1;
        indices{2} = n2;
        r=PhasorArray(pA1(indices{:}),"reduce", false);
        if nargin == 1
            r = reshape(r,size(n1,1),size(n1,2));
        end
    end

    function r = vect(pA1)
        %VECT transform phasorArray matrix to column vector (col operator)
        %   This function takes a PhasorArray object and transforms it into
        %   a column vector by stacking the columns of the input matrix.
        %
        %   Syntax:
        %       r = vect(pA1)
        %
        %   Input:
        %       pA1 - PhasorArray object to be transformed.
        %
        %   Output:
        %       r  - Column vector (PhasorArray) obtained by stacking the
        %            columns of the input matrix.
        %
        %   Example:
        %       % Assuming pA1 is a PhasorArray object
        %       r = vect(pA1);
        %
        %   See also: reshape, pvalue
        pval = pvalue(pA1);
        dims = size(pval);
        dimsCell = num2cell(dims);
        r1 = reshape(pval, [], 1, dimsCell{3:end});
        r = PhasorArray(r1, "reduce", false);
    end
    function r = pad(pA1,delta_h)
        %PAD the PhasorArray with delta_h phasor
        %   - r = pad(pA1,delta_h) : pad the PhasorArray pA1 with delta_h 0 phasor
        %   - r = pad(pA1,[delta_h1 delta_h2 delta_h3]) : pad the PhasorArray pA1 with delta_h1 0 phasor along the first dimension, delta_h2 0 phasor along the second dimension and delta_h3 0 phasor along the third dimension
        %
        %pad (add zeros) to phasor array
        %pad(A, h) add h phasor to the phasor array
        %pad(A, [h1 h2 h3]) pad A in each direction
        r=PhasorArrayPad(pA1,delta_h);
    end
    function r = padh(pA1,delta_h)
        %PAD the PhasorArray with delta_h phasor
        %   - r = pad(pA1,delta_h) : pad the PhasorArray pA1 with delta_h 0 phasor
        o1v = pA1.value;
        o1v = cat(3,zeros(size(o1v,1),size(o1v,2),delta_h),o1v,zeros(size(o1v,1),size(o1v,2),delta_h));
        r=PhasorArray(o1v);
    end

    function r = ctranspose(pA1)
        %CTRANSPOSE overloading for PhasorArray
        %CTRANSPOSE(A) : A' = A(t)'
        r = mctranspose(pA1);
    end
    function r = transpose(pA1)
        %TRANSPOSE overloading for PhasorArray
        %TRANSPOSE(A) : A.' = A(t).'
        r = mtranspose(pA1);
    end
    function r = pagectranspose(pA1)
        %PAGECTRANSPOSE overloading for PhasorArray
        %PAGECTRANSPOSE(A) apply ctranspose to each page of A resulting in A(-t)'
        r = PhasorArray(pagectranspose(pvalue(pA1)));
    end
    function r = pagetranspose(pA1)
        %PAGETRANSPOSE overloading for PhasorArray
        %PAGETRANSPOSE(A) apply transpose to each page of A resulting in A(t).'
        r = PhasorArray(pagetranspose(pvalue(pA1)));
    end
    function r = mtranspose(pA1)
        % mtranspose(A(t))
        r=pagetranspose(pA1);
    end
    function r = mctranspose(pA1)
        % mctranspose(A(t))
        r=PhasorArray(flip(pagectranspose(pA1.value),3));
    end

    function r = pmax(pA1,pA2,nvp)
        % pmax Compute elementwise or phasorwise max between two PhasorArray objects.
        %
        %   r = pmax(pA1, pA2, 'method', method) computes the maximum between two
        %   PhasorArray objects `pA1` and `pA2` according to the specified method.
        %
        %   Inputs:
        %       pA1, pA2  - PhasorArray objects to compare.
        %       method   - (Optional) Specifies the max computation method:
        %                  * 'elementwise' (default): Max computed based on norm.
        %                  * 'phasorwise': Max computed independently for each phasor.
        %
        %   Output:
        %       r       - The resulting PhasorArray after applying the max operation.
        %
        %   Notes:
        %   - In 'elementwise' mode, each phasor's norm is computed first, and the
        %     element-wise max is applied to select the dominant phasor.
        %   - In 'phasorwise' mode, each corresponding phasor component is compared directly.
        %   - Both inputs are unified in size using `PhasorUnif` before processing.
        %
        %   See also: PhasorUnif, pvalue, norm
        arguments
            pA1
            pA2
            nvp.method char {mustBeMember(nvp.method,{'elementwise','phasorwise'})} = 'elementwise'
        end
        if strcmp(nvp.method,'elementwise')
            d=PhasorUnif(pA1,pA2);
            pA1 = d{1};
            pA2 = d{2};

            % norm of each phasor (along the third dim)
            n1=sum(pA1.value .* conj(pA1.value),3);
            n2=sum(pA2.value .* conj(pA2.value),3);
            %max of the norm
            n=max(n1,n2);


            %elementwise max
            r=PhasorArray(pvalue(pA1).*(n1==n) + pvalue(pA2).*(n1~=n));
        else
            d=PhasorUnif(pA1,pA2);
            r = PhasorArray(max(pvalue(d{1}),pvalue(d{2})));
        end
    end
    
    function out = dc(pA1)
        %DC Constant part: harmonic 0 alone, i.e. the mean of A(t) over one period.
        %
        %   See also: ac, phas.
        out = PhasorArray(pA1{:,:,0});
    end


    function out = ac(pA1)
        %AC Oscillating part: every harmonic except 0, so A - dc(A).
        %
        %   See also: dc, phas.
        negp = pA1{:,:,-pA1.h:-1};
        posp = pA1{:,:,1:pA1.h};
        out = PhasorArray(cat(3,negp,zeros(size(pA1,1),size(pA1,2)),posp));
    end

    function varargout=phas(pA1,h)
        % phas Extract the phasor of order `h` from a PhasorArray.
        %
        %   phas(pA1, h) returns the phasor of order `h` from the PhasorArray `pA1`.
        %
        %   Inputs:
        %       pA1 - The PhasorArray object.
        %       h  - The order of the phasor to extract.
        %            - If `h` is a scalar, returns a single phasor slice.
        %            - If `h` is a vector, returns multiple slices at the specified orders.
        %
        %   Outputs:
        %       - A 3D array corresponding to the phasor(s) of order `h`.
        %
        %   Behavior:
        %   - If `abs(h) <= pA1.h`, extracts the corresponding phasor.
        %   - If `h` is **out of range**, pads `pA1` to match the highest requested order.
        %   - If `h` is **beyond padding limits**, returns a zero matrix.
        %
        %   See also: pad, reduce, PhasorArray
        if strcmp(h,':')
            varargout{1} = pA1.value;
            return
        end
        if isscalar(h)
            if abs(h)<=pA1.h
                indices = repmat({':'}, 1, ndims(pA1));
                indices{3} = (size(pA1,3)+1)/2+h;
                varargout{1} = pA1(indices{:});
            else
                sz = size(pA1);
                sz(3) = 1;
                varargout{1} = zeros(sz);
            end
        else
            if max(abs(h))>=pA1.h
                pA1=pA1.pad( max(abs(h))-pA1.h);
            end
            varargout{1}=pA1(:,:,(end+1)/2+h);
        end

    end
    function obj =phasAssign(obj,h,varargin)
        % phasAssign Assign values to specific phasor components in a PhasorArray.
        %
        %   obj = phasAssign(obj, h, values) returns a copy of the PhasorArray
        %   with the phasor at index h replaced by the provided values.
        %
        %   Inputs:
        %       obj    - The PhasorArray object to modify.
        %       h      - The index (or indices) of the phasor(s) to be assigned.
        %       values - The new values to be assigned to the specified phasors.
        %
        %   Outputs:
        %       obj    - A modified PhasorArray with updated phasor values.
        %
        %   Notes:
        %   - If h exceeds the current order of the PhasorArray, it is automatically expanded.
        %   - This method does not modify the original object but returns a modified copy.
        %   - The provided values must be compatible in size with the existing array dimensions.
        %
        %   See also: phas, sub, reduce, trunc
        obj{:,:,h}=varargin{:};
    end


    % --- Comparisons act on A(t), not on the coefficient array: entry (i,j)
    %   judges A_ij(t) against B_ij(t) over a period. Comparing coefficients
    %   would fail silently, MATLAB's < testing real parts only on complex.
    function [r, frac, crenel, Cph] = lt(pA1,pA2)
        %LT  A < B, meaning A_ij(t) < B_ij(t) for every t. Returns [n x m].
        %   [R, FRAC, CRENEL, CPH] = LT(A,B) grades the answer: FRAC is the
        %   fraction of the period over which it holds, CRENEL the indicator
        %   over theta, CPH the same as a PhasorArray. A < B takes R alone.
        %   See also gt, le, ge, eq, compareInTime.
        [r, frac, crenel, Cph] = compareInTime(pA1, pA2, @lt);
    end
    function [r, frac, crenel, Cph] = gt(pA1,pA2)
        %GT  A > B, meaning A_ij(t) > B_ij(t) for every t. Returns [n x m].
        %   Extra outputs as in LT.
        %   See also lt, le, ge, eq, compareInTime.
        [r, frac, crenel, Cph] = compareInTime(pA1, pA2, @gt);
    end
    function [r, frac, crenel, Cph] = le(pA1,pA2)
        %LE  A <= B, meaning A_ij(t) <= B_ij(t) for every t. Returns [n x m].
        %   Extra outputs as in LT.
        %   See also lt, gt, ge, eq, compareInTime.
        [r, frac, crenel, Cph] = compareInTime(pA1, pA2, @le);
    end
    function [r, frac, crenel, Cph] = ge(pA1,pA2)
        %GE  A >= B, meaning A_ij(t) >= B_ij(t) for every t. Returns [n x m].
        %   Extra outputs as in LT.
        %   See also lt, gt, le, eq, compareInTime.
        [r, frac, crenel, Cph] = compareInTime(pA1, pA2, @ge);
    end
    function r = abs(pA1, nvp)
        %ABS  |A(t)| taken pointwise in time, as a PhasorArray.
        %   ABS(A) samples A over one period, takes the modulus and transforms
        %   back, returning everything the sampling supports. ABS(A, h=H)
        %   truncates to order H; TRUNC does the same afterwards.
        %
        %   |A(t)| has a corner wherever A crosses zero, so it is not band
        %   limited. Being continuous it does not ring, but it converges only
        %   as 1/h: measured relative error 41% at h = A.h, 8% at 4*A.h and
        %   1.2% at the full sampling order. Truncating back to A.h is
        %   therefore a real approximation, not a formality.
        %
        %   Works on complex A(t), where the modulus is still real.
        %
        %   See also lt, evalp, trunc.
        arguments
            pA1
            nvp.h {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
        end
        % Sampled then transformed back, so it needs numbers: |A(t)| is not an
        % expression a decision variable can carry, and its Fourier series is
        % not the one of A anyway -- the modulus is not analytic.
        solvedPayload(pvalue(pA1), "abs");
        Nt = 2^nextpow2(max(256, 16*(pA1.h+1)));
        th = (0:Nt-1)/Nt * 2*pi;
        r  = PhasorArray(TimeArray2Phasors(abs(evalp(pA1, th)), 1, th));
        if ~isempty(nvp.h), r = r.trunc(nvp.h); end
    end
    function r = regularize(pA1, epsilon, nvp)
        %REGULARIZE  Mollify A(t) so its Fourier series converges uniformly.
        %   REGULARIZE(A, EPS) returns A_eps = phi_eps * A, the convolution
        %   with the C^inf bump phi(t) = exp(-1/(1-t^2)) on |t|<1, scaled to
        %   width EPS and normalised to unit mass. EPS is a time, in the same
        %   units as T (name-value, default 2*pi).
        %
        %   This is the regularization of the erratum to "A TBLMI Framework
        %   for Harmonic Robust Control" (Vernerey, Riedinger, Daafouz).
        %   Banded truncation converges in operator norm only under uniform
        %   convergence of the Fourier series, which A in L^inf alone does not
        %   give; mollifying yields a C^inf A_eps for which it does, and the
        %   trajectories converge uniformly back as EPS -> 0.
        %
        %   Convolution in time is multiplication of the coefficients, so this
        %   costs one windowing: A_k becomes phihat(EPS*k*omega)*A_k. phi is
        %   real and even, so the window is real -- no phase is introduced.
        %
        %   Exact, with no edge effect to worry about: A is periodic, hence
        %   defined on all of R, so there is no finite record to convolve over
        %   and nothing to extend. Checked to 1e-15 against an explicit
        %   convolution sampled over [-2*pi, 4*pi], for EPS up to 4.
        %
        %   Apply it where the result is knowingly non-smooth, not by default
        %   after every transform: on band-limited data it removes nothing and
        %   costs O(EPS^2) of bias, and it does not undo aliasing -- a fold
        %   deposits energy at low k, where the window is 1. Oversampling is
        %   the defence against that, and abs and the comparisons already do it.
        %
        %   phihat decreases from 1 to 0 over [0, 5] and first vanishes at
        %   4.9965, its envelope then falling sub-exponentially: 4.6e-02 near
        %   10, 6.8e-04 near 40. Keeping EPS*h*omega inside that first lobe
        %   leaves the window positive and monotone; past it the window has
        %   zeros, which is faithful mollification but no longer reads as a
        %   low-pass.
        %
        %   See also abs, trunc, neglect.
        arguments
            pA1
            epsilon (1,1) double {mustBeNonnegative}
            nvp.T   (1,1) double {mustBePositive} = 2*pi
        end
        if epsilon == 0, r = pA1; return, end
        k = reshape(-pA1.h:pA1.h, 1, 1, []);
        % The weight is one number per harmonic, grown to the payload's size:
        % sym and sdpvar do not broadcast a 1x1xK against an n x m x K. It is
        % not a phasor product either -- it scales coefficient k, it does not
        % convolve with it.
        w = mollifierFT(epsilon * k * (2*pi/nvp.T));
        v = pvalue(pA1);
        r = PhasorArray(v .* repmat(w, size(v, 1), size(v, 2), 1));
    end
    function r = eq(pA1,pA2)
        %EQ  A == B, meaning A_ij(t) = B_ij(t) for every t. Returns [n x m].
        %
        %   No sampling is involved: the map between a periodic signal and its
        %   harmonic coefficients is a bijection, so equality of the functions
        %   is exactly equality of the coefficients once MINUS has padded both
        %   to a common order. For a tolerance, use iszero(A-B, tol).
        %
        %   See also ne, iszero, lt.
        % Only decision variables are refused: == 0 on one builds a constraint
        % rather than a logical, so all(...) would have nothing to reduce, and
        % equality of unsolved variables is a constraint to state, not a
        % question to ask. sym answers this perfectly well and is left alone.
        dv = pvalue(pA1 - pA2);
        if isa(dv, 'ndsdpvar') || isa(dv, 'sdpvar')
            error('PhasorArray:eq:decisionVariable', ...
                ['== compares values and this array holds unsolved decision ' ...
                 'variables. Solve first and compare sdpval of the result, or ' ...
                 'state the equality as a constraint.']);
        end
        r = all(dv == 0, 3);
    end
    function r = ne(pA1,pA2)
        %NE  A ~= B, the negation of EQ. Returns [n x m].
        %   See also eq.
        r = ~eq(pA1,pA2);
    end
    function r = double(pA1)
        %DOUBLE Convert PhasorArray to double precision.
        %   r = DOUBLE(pA1) converts the PhasorArray object pA1 to a double precision array.
        r= double(pvalue(pA1));
    end


    function pA2 = repeat(pA1,m)
        %REPEAT Repeat the PhasorArray object over its period by a factor m
        %
        %   pA2 = repeat(pA1, m) repeats the PhasorArray object pA1 over its
        %   period by a factor m. If m is negative, the signal is reversed
        %   in time, effectively dividing the period by m.
        %
        %   Inputs:
        %       pA1 - The PhasorArray object to be repeated.
        %       m  - An integer factor by which to repeat the period.
        %            Default value is 2. If m is negative, the signal is
        %            reversed in time.
        %
        %   Outputs:
        %       pA2 - The resulting PhasorArray object after repeating the
        %            period by a factor m.
        %
        %   Example:
        %       pA1 = PhasorArray(...); % Create a PhasorArray object
        %       pA2 = repeat(pA1, 3);    % Repeat the period of pA1 by a factor of 3
        %
        %   See also: PhasorArray.zeros, retro
        arguments
            pA1
            m {mustBeInteger(m)} = 2
        end
        h=pA1.h;
        h_new=h*abs(m);
        pA2 = PhasorArray.zeros(size(pA1,1),size(pA1,2),2*h_new+1);
        pA2(:,:,1:abs(m):end)=pA1.value;
        if m<0
            pA2=retro(pA2);
        end

    end


    function r = antiD(pA1,T)
        %ANTI-D Compute the anti-derivative (primitive) of non-zero phasors of a PhasorArray.
        %   r = ANTI-D(pA1, T) computes the anti-derivative of the PhasorArray pA1
        %   with respect to the period T. The function ignores the 0-th phasor to
        %   produce a periodic output.
        %
        %   INPUTS:
        %       pA1 - The PhasorArray object to be anti-differentiated.
        %       T  - (Optional) The period of the PhasorArray. Defaults to 2*pi.
        %
        %   OUTPUT:
        %       r  - The resulting PhasorArray after computing the anti-derivative.
        %
        %   EXAMPLE:
        %       A = PhasorArray.random(3, 3, 5);
        %       B = A.antiD(2*pi);
        %
        %   NOTE:
        %       The function divides each non-zero phasor by (1i * k * 2 * pi / T),
        %       where k is the phasor index. The 0-th phasor is set to zero.
        arguments
            pA1
            T=2*pi
        end
        pA1=pA1.value;
        [~,~,n3]=size(pA1);
        for ii=1:n3
            k=-(n3-1)/2+ii-1;
            if k~=0
                pA1(:,:,ii)=pA1(:,:,ii)/(1i*k*2*pi/T);
            else
                pA1(:,:,ii)=0;
            end
        end
        r=PhasorArray(pA1);
    end

    function r = d(pA1,T)
        % d - Derive phasor array with respect to a given period
        %
        % Syntax: r = d(pA1, T)
        %
        % Inputs:
        %    pA1 - The input phasor array object
        %    T  - The period with respect to which the phasor array is derived (default is 2*pi)
        %
        % Outputs:
        %    r  - The resulting phasor array after derivation
        %
        % Example:
        %    phasorArray = PhasorArray(someValues);
        %    result = phasorArray.d(phasorArray, 2*pi);
        arguments
            pA1
            T=2*pi
        end
        %[~,~,n3]=size(pA1);
        try
            K=permute((-pA1.h):(pA1.h),[1 3 2])*1i*2*pi/T;
            oo=bsxfun(@times,pA1.value,K);
            r=PhasorArray(oo);
        catch e
            % warning(e.message)
            % warning('Attempting to diff a non basic type, switching to manual diff')
            % We need to manually perform the operation intended for each "slice" of P.
            for idx = 1:size(pA1, 3)
                pA1(:,:,idx) = pA1(:,:,idx) * K(idx);
            end
            r = pA1;
            % warning('done with success')
        end
    end
    function r = PhaseShift(pA1,angle)
        % PhaseShift Apply a phase shift to a PhasorArray.
        %
        %   r = PhaseShift(pA1, angle) shifts the phase of the PhasorArray `pA1`
        %   by `angle` radians.
        %
        %   Inputs:
        %       pA1    - The PhasorArray object.
        %       angle - The phase shift value(s). It can be:
        %               - A scalar (applies to all elements of `pA1`).
        %               - A row vector (broadcasts across `pA1` columns if `pA1` is a column).
        %               - A column vector (broadcasts across `pA1` rows if `pA1` is a row).
        %               - A matrix (applies element-wise if `pA1` is scalar).
        %
        %   Output:
        %       r - The phase-shifted PhasorArray.
        %
        %   Behavior:
        %   - If `pA1` is a scalar, `angle` can be a **matrix**, applying element-wise shifts.
        %   - If `pA1` is a **row vector**, `angle` must be a column vector (broadcasted).
        %   - If `pA1` is a **column vector**, `angle` must be a row vector (broadcasted).
        %   - If `pA1` is a **matrix**, `angle` must be scalar.
        %
        %   Errors:
        %   - If `pA1` and `angle` are both row or column vectors, the function raises an error.
        %   - If `pA1` is a matrix, `angle` must be a scalar.
        %
        %   See also: dephase, PhasorArray
        if numel(angle)>1
            if isrow(pA1) && iscolumn(angle)
                n1=numel(angle);
                %n2=pA1.size(2);
                r=repmat(pA1,n1,1);
                for angli=1:numel(angle)
                    angle(angli);
                    r{angli,:}=pA1.PhaseShift(angle(angli));
                end
                return
            elseif iscolumn(pA1) && isrow(angle)
                %n1=pA1.size(1);
                n2=numel(angle);
                r=repmat(pA1,1,n2);
                for angli=1:numel(angle)
                    r{:,angli}=pA1.PhaseShift(angle(angli));
                end
                return
            elseif isscalar(pA1) && ismatrix(angle)
                [n1,n2]=size(angle,[1,2]);

                r=repmat(pA1,n1,n2);
                for angli=1:size(angle,1)
                    for anglj=1:size(angle,2)
                        r{angli,anglj}=pA1.PhaseShift(angle(angli,anglj));
                    end
                end
                return
            elseif iscolumn(pA1) && iscolumn(angle)
                error('PhasorArray:phaseShift:invalidAngleShape', 'If pA1 is a column vector, angle must be a row vector (and vice versa).')
            elseif isrow(pA1) && isrow(angle)
                error('PhasorArray:phaseShift:invalidAngleShape', 'If pA1 is a row vector, angle must be a column vector (and vice versa).')
            else
                error('PhasorArray:phaseShift:invalidAngleShape', 'angle must be a scalar for a matrix-valued PhasorArray.')
            end
        else
            r=dephase(pA1,angle);
        end
    end
    function [oInv, oInvt, norm_err, norm_ref] = inv(pA1, nvp)
        %INV Compute the phasor representation of the pointwise inverse of A(t).
        %
        %   This function computes the **pointwise inverse** of the time-domain realization
        %   of the periodic matrix A(t) and then reconstructs its phasor representation.
        %
        %   Syntax:
        %   r = INV(pA1)
        %       Computes the phasors of the pointwise inverse of `pA1` using default parameters.
        %
        %   Name-Value Arguments (forwarded to PhasorInv):
        %       'nT'              - Number of periods for time-domain evaluation (default: 1).
        %       'T'               - Period used for simulation (default: 1).
        %       'm'               - Log2 of time-domain discretization points.
        %       'plot'            - Plot A⁻¹(t) after computation (default: false).
        %       'reduceThreshold' - Threshold for phasor reduction (default: 4e-15).
        %       'verbose'         - Print debug information (default: false).
        %
        %   See also: PhasorInv, mldivide, mrdivide
        arguments
            pA1  PhasorArray
            nvp.nT              = 1
            nvp.T               = 2*pi
            nvp.m               = []
            nvp.plot            = false
            nvp.reduceThreshold = 4e-15
            nvp.reduceMethod    = 'relative'
            nvp.autoTrunc       = false
            nvp.verbose         = false
            nvp.evalInv         = false
        end

        C = namedargs2cell(nvp);

        % Propagate nargout explicitly: PhasorInv conditionally computes
        % norm_err / norm_ref only when nargout > 2 (expensive computation).
        if nargout <= 1
            oInv = PhasorArray(PhasorInv(pA1, C{:}));
        elseif nargout == 2
            [Ainvph, oInvt] = PhasorInv(pA1, C{:});
            oInv = PhasorArray(Ainvph);
        else
            [Ainvph, oInvt, norm_err, norm_ref] = PhasorInv(pA1, C{:});
            oInv = PhasorArray(Ainvph);
        end
    end
    function [PhDet,det_t] = det(pA1,nvp)
        %DET Compute the pointwise determinant of A(t) and reconstruct its phasors.
        %
        %   This function computes the determinant of the **time-domain realization** A(t),
        %   then reconstructs its phasor representation via Fourier Transform.
        %
        %   The determinant is computed **pointwise at each time step**, meaning:
        %   - If A(t) is periodic, det(A(t)) is computed over time.
        %   - The final result represents the phasor decomposition of det(A(t)).
        %
        %   Syntax:
        %   [PhDet, det_t] = DET(pA1)
        %       Computes the phasors of det(A(t)) using default settings.
        %
        %   [PhDet, det_t] = DET(pA1, 'nT', nT, 'T', T, 'm', m,
        %                             'plot', plotFlag, 'autoTrunc', autoTrunc)
        %       Computes det(A(t)) with additional control over truncation and plotting.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The PhasorArray object whose determinant is computed.
        %
        %   Name-Value Pair Arguments:
        %   - 'nT' (integer, optional) : Number of periods used in the time-domain evaluation. Default: 1.
        %   - 'T' (double, optional) : The period used for simulation. Default: 1.
        %   - 'm' (integer, optional) :
        %       - Power of two controlling time-domain discretization.
        %       - Can be set to [] for automatic selection based on the number of phasors.
        %   - 'plot' (logical, optional) : If true, plots det(A(t)) after computation. Default: false.
        %   - 'autoTrunc' (logical, optional) :
        %       - true : Uses the derivative of phasors to **automatically detect**
        %         the significant number of phasors.
        %       - false (default) : Uses a fixed threshold-based reduction method.
        %
        %   - If 'autoTrunc' is false, the following options apply:
        %       - 'reduceThreshold' (double, optional) : The threshold for reducing phasors. Default: 1e-20.
        %       - 'reduceMethod' (char, optional) : Reduction strategy.
        %           - 'absolute' : Remove phasors with magnitude < reduceThreshold.
        %           - 'relative' (default) : Remove phasors with magnitude < max(magnitude) * reduceThreshold.
        %
        %   Output Arguments:
        %   - PhDet (PhasorArray) : The phasor representation of det(A(t)).
        %   - det_t (array) : The time-domain realization of det(A(t)).
        %
        %   Algorithm:
        %   1. Compute A(t) in the time domain using an **IFFT**.
        %   2. Perform **pointwise determinant computation** det(A(t)).
        %   3. Reconstruct phasors by applying an **FFT** to det(A(t)).
        %
        %   Example:
        %   % Compute the phasors of det(A(t)) using default settings
        %   [PhDet, det_t] = det(A);
        %
        %   % Compute det(A(t)) over 2 periods with auto truncation
        %   [PhDet, det_t] = det(A, 'nT', 2, 'autoTrunc', true);
        %
        %   % Compute det(A(t)) with manual truncation and thresholding
        %   [PhDet, det_t] = det(A, 'reduceThreshold', 1e-15, 'reduceMethod', 'absolute');
        %
        %   See also: INV, REDUCE, FFT, IFFT.
        arguments
            pA1
            nvp.nT=1
            nvp.T=2*pi
            nvp.m=[]
            nvp.plot=false
            nvp.reduceThreshold = 1e-20
            nvp.reduceMethod = 'relative'
            nvp.autoTrunc = false
        end


        C=namedargs2cell(nvp);
        [PhDet,det_t] =  PhasorDet(pA1,C{:});
        PhDet=PhasorArray(PhDet);

    end

    function r = diag(pA1,K)
        %DIAG Extract or construct a diagonal PhasorArray.
        %
        %   This function operates similarly to MATLAB's `diag`:
        %   - If `pA1` is a **vector**, it constructs a **diagonal PhasorArray**.
        %   - If `pA1` is a **square matrix**, it extracts its diagonal as a **vector**.
        %
        %   Syntax:
        %   r = DIAG(A)
        %       Constructs a diagonal PhasorArray from a vector `A`.
        %
        %   r = DIAG(A, K)
        %       Constructs a diagonal PhasorArray, with `K` as the diagonal offset.
        %
        %   r = DIAG(A)
        %       Extracts the diagonal as a PhasorArray vector when `A` is square.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : Input PhasorArray, either a vector or a square matrix.
        %   - K (integer, optional) : Offset for the diagonal placement. Default: 0.
        %
        %   Output:
        %   - r (PhasorArray) :
        %       - If `pA1` is a vector, `r` is a **diagonal PhasorArray**.
        %       - If `pA1` is square, `r` is the **extracted diagonal**.
        %
        %   Example:
        %   % Create a diagonal PhasorArray from a vector
        %   A = PhasorArray(rand(5,1));
        %   R = diag(A);
        %
        %   % Extract the diagonal from a square PhasorArray
        %   B = PhasorArray(rand(4,4));
        %   d = diag(B);
        %
        %   % Construct a diagonal PhasorArray with an offset
        %   C = diag(A, 1);
        %
        %   See also: TRACE.
        arguments
            pA1
            K=0
        end
        if isscalar(pA1) && K==0
            r = pA1;
            return
        end

        if isvector(pA1)
            if ~iscolumn(pA1)
                pA1=pA1.';
            end
            if isa(pA1.value,"sdpvar")
                r = PhasorArray(diag(pA1));
                return
            end

            if isspecial(pA1)
                for ii = 1:size(pA1,3)
                    r{ii}  = PhasorArray(diag(pA1(:,:,ii)));
                end
                r=cat(3,r{:});
                return
            end

            r = PhasorArray.zeros(numelt(pA1)+abs(K));


            I = logical(diag(ones(numelt(pA1),1),K));
            r{I}=r{I}+pA1;
            return
        end
        assert(issquare(pA1),"diag function : PhasorArray must represent a square matrix");
        n=size(pA1,1);
        I=(1:n)*(n+1)-n;
        r=pA1{I};
    end
    function r = trace(pA1,type)
        %TRACE Compute the phasor representation of the trace of A(t).
        %
        %   This function computes the **trace** of a PhasorArray A, which is the sum
        %   of its diagonal elements. It supports two modes:
        %   - 'phasor' (default) : Returns the **full** phasor representation of the trace.
        %   - '0'                : Returns only the **DC (0th order)** component of the trace.
        %
        %   Syntax:
        %   r = TRACE(A)
        %       Computes the phasor representation of the trace of `A`.
        %
        %   r = TRACE(A, 'phasor')
        %       Same as above; returns the full phasor trace.
        %
        %   r = TRACE(A, '0')
        %       Returns only the DC component of the trace.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : A **square PhasorArray**.
        %   - type (char, optional) :
        %       - 'phasor' (default) : Returns the **full phasor representation** of trace(A).
        %       - '0' : Returns only the **DC component** (0th phasor) of trace(A).
        %
        %   Output:
        %   - r (PhasorArray or numeric) :
        %       - If 'phasor', `r` is a **PhasorArray** containing the full phasor trace.
        %       - If '0', `r` is a **scalar numeric value** (the DC component of trace(A)).
        %
        %   Example:
        %   % Compute the full phasor trace of a PhasorArray
        %   A = PhasorArray(rand(3,3,11));
        %   tr_phasor = trace(A);
        %
        %   % Compute only the DC component of the trace
        %   tr_DC = trace(A, '0');
        %
        %   See also: DIAG.
        arguments
            pA1
            type char {mustBeMember(type,{'phasor','0'})} = 'phasor'
        end
        r=sum(diag(pA1));
        if strcmp(type,'0')
            %return the DC value of the trace
            r=r{:,:,0};
        end
    end

    function r = numelt(pA1)
        %NUMELT Compute the number of elements in the first two dimensions of A.
        %
        %   This function returns `size(A,1) * size(A,2)`, which represents
        %   the number of matrix elements per time instant.
        %
        %   Syntax:
        %   r = NUMELT(A)
        %       Computes the total number of elements in A(t).
        %
        %   Input:
        %   - pA1 (PhasorArray) : A PhasorArray object.
        %
        %   Output:
        %   - r (integer) : Number of elements in `A(t)`.
        %
        %   Example:
        %   % Compute the number of elements in a 4x5 PhasorArray
        %   A = PhasorArray(rand(4,5,10));
        %   e = numelt(A);  % Returns 4*5 = 20.
        r=size(pA1,1)*size(pA1,2);
    end
    function r = h(pA1)
        %H Compute the maximal phasor order stored in A.
        %
        %   This function extracts the **highest phasor order** present in a PhasorArray.
        %
        %   Syntax:
        %   r = H(A)
        %       Computes the highest stored harmonic order.
        %
        %   Input:
        %   - pA1 (PhasorArray) : A PhasorArray object.
        %
        %   Output:
        %   - r (integer) : The highest stored harmonic order in A.
        %
        %   Example:
        %   % Get the maximum harmonic order in a PhasorArray
        %   A = PhasorArray(rand(4,4,11));
        %   max_h = h(A);  % Returns (11-1)/2 = 5.
        %
        %   See also: DIM.
        r=(size(pA1,3)-1)/2;
    end
    function r = dim(pA1)
        %DIM Compute the number of elements in A(t).
        %
        %   This function returns `size(A,1) * size(A,2)`, which represents
        %   the number of matrix elements per time instant.
        %
        %   Syntax:
        %   r = DIM(A)
        %       Computes the total number of elements in A(t).
        %
        %   Input:
        %   - pA1 (PhasorArray) : A PhasorArray object.
        %
        %   Output:
        %   - r (integer) : Number of elements in `A(t)`.
        %
        %   Example:
        %   % Compute the number of elements in a 4x5 PhasorArray
        %   A = PhasorArray(rand(4,5,10));
        %   e = dim(A);  % Returns 4*5 = 20.
        %
        %   See also: NUMELT.
        r=size(pA1,1)*size(pA1,2);
    end

    function r = T_bt(pA1, m, h2)
        %T_BT Block-Toeplitz operator of a PhasorArray, harmonics ascending.
        %
        %   T_bt is the canonical name of this operator, matching the notation
        %   used in the publications and the naming of N_bt and S_bt. BT is the
        %   original name, kept as an alias.
        %
        %   Syntax:
        %   r = T_bt(A, h)        - Square, size (2*h+1)N x (2*h+1)N.
        %   r = T_bt(A, h1, h2)   - Rectangular, size (2*h1+1)N x (2*h2+1)N.
        %
        %   See also: BT, T_tb, N_bt, S_bt.
        arguments
            pA1
            m = []
            h2 = []
        end
        r = BT(pA1, m, h2);
    end

    function r = T_tb(pA1, m, h2)
        %T_TB Toeplitz-blocks operator of a PhasorArray, harmonics ascending.
        %
        %   T_tb is the canonical name of this operator, matching the notation
        %   used in the publications and the naming of N_tb and S_tb. TB is the
        %   original name, kept as an alias.
        %
        %   Syntax:
        %   r = T_tb(A, h)        - Square, size (2*h+1)N x (2*h+1)N.
        %   r = T_tb(A, h1, h2)   - Rectangular, size (2*h1+1)N x (2*h2+1)N.
        %
        %   Note that T_tb(A*B) differs from T_tb(A)*T_tb(B) at finite h: the
        %   product loses the harmonics that truncation discards.
        %
        %   See also: TB, T_bt, N_tb, S_tb.
        arguments
            pA1
            m = []
            h2 = []
        end
        r = TB(pA1, m, h2);
    end

    function r = BT(pA1, m, h2)
        %BT Construct a Block Toeplitz matrix from a PhasorArray.
        %
        %   This function builds a **Block Toeplitz matrix** of size (2*h1+1)N x (2*h2+1)N
        %   from the given PhasorArray `A`.
        %
        %   Syntax:
        %   r = BT(A, h)         - Square block matrix of size (2*h+1)N x (2*h+1)N.
        %   r = BT(A, h1, h2)    - Rectangular block matrix of size (2*h1+1)N x (2*h2+1)N.
        %   r = BT(A, [h1, h2])   - Rectangular block matrix using a vector input.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The input PhasorArray.
        %   - m (integer, optional) : Harmonic order h (or vector [h1, h2]).
        %   - h2 (integer, optional) : Second harmonic order for rectangular truncation.
        %
        %   See also: TB, spBT, array2BToeplitz.
        arguments
            pA1
            m = []
            h2 = []
        end
        if nargin > 2
            m = [m, h2];
        end
        r = array2BToeplitz(pA1, m);
    end
    function r = TB(pA1, m, h2)
        %TB Construct a Toeplitz Block matrix from a PhasorArray.
        %
        %   This function constructs a **Toeplitz Block (TB) matrix** of size (2*h1+1)N x (2*h2+1)N
        %   from the given PhasorArray `A`.
        %
        %   Syntax:
        %   r = TB(A, h)         - Square block matrix of size (2*h+1)N x (2*h+1)N.
        %   r = TB(A, h1, h2)    - Rectangular block matrix of size (2*h1+1)N x (2*h2+1)N.
        %   r = TB(A, [h1, h2])   - Rectangular block matrix using a vector input.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The input PhasorArray.
        %   - m (integer, optional) : Harmonic order h (or vector [h1, h2]).
        %   - h2 (integer, optional) : Second harmonic order for rectangular truncation.
        %
        %   See also: BT, spTB, array2TBlocks.
        arguments
            pA1
            m = []
            h2 = []
        end
        if nargin > 2
            m = [m, h2];
        end
        r = array2TBlocks(pA1, m);
    end
    function r = spTB(pA1, m, h2)
        %SPTB Construct a sparse Toeplitz Block representation of A.
        %
        %   Syntax:
        %   r = spTB(A, h1, h2)
        %
        %   See also: TB, BT, sparray2TBlocks.
        arguments
            pA1
            m = []
            h2 = []
        end
        if nargin > 2
            m = [m, h2];
        end
        r = sparray2TBlocks(pA1, m);
    end
    function r = spBT(pA1, m, h2)
        %SPBT Sparse Block-Toeplitz operator of a PhasorArray.
        %
        %   Same result as BT, returned as a sparse matrix. The operator is
        %   banded whenever the array is truncated, so the sparse form is worth
        %   using as soon as the harmonic order grows.
        %
        %   See also: BT, T_bt, spTB.
        arguments
            pA1
            m = []
            h2 = []
        end
        if nargin > 2
            m = [m, h2];
        end
        r = sparse(pA1.BT(m));
    end
    function r = FvTB(pA1,h)
        %FVTB Compute the Fourier representation of the vectorized form of A(t).
        %
        %   This function first **vectorizes** the time-dependent matrix A(t) by stacking
        %   all its columns into a **single-column vector** a(t) = vect(A(t)), and then
        %   applies the **F_tb** transformation to obtain its Fourier representation.
        %
        %   **Alternative Computation:**
        %   - In MATLAB notation, vectorizing A(t) is equivalent to **A{:}**.
        %   - Thus, this function is **equivalent to applying F_tb to A{:}**, i.e.:
        %       FvTB(A, h) = F_tb(A{:}, h).
        %
        %   **Procedure:**
        %   1. Compute the column-wise vectorization:
        %       - If A(t) is an `N×M` matrix, then vect(A(t)) is a column vector of size `NM × 1`.
        %       - This stacking is done **column by column** (following MATLAB’s column-major order).
        %   2. Compute the Fourier series coefficients **up to order h**:
        %       - Instead of applying `F_tb` directly to A(t), we apply it to a(t) = vect(A(t)).
        %       - The result is the Fourier representation of `vect(A(t))`, which is a vectorized form
        %         of `F_tb(A, h)`.
        %
        %   **Key Property:**
        %   If `y(t) = A(t)x(t)`, then:
        %       FvTB(y) = FvTB(A) ⋅ FvTB(x),
        %   where `FvTB(A, h)` is the **Fourier vectorized representation** of A.
        %
        %   **Dimension of the Output:**
        %   - If `A` is an `N×M` matrix, then `FvTB(A, h)` is a `((2h+1)NM × 1)` vector.
        %
        %   Syntax:
        %   r = FVTB(A, h)
        %       Computes the Fourier representation of `vect(A)` up to order `h`.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The input PhasorArray representing A(t).
        %   - h (integer) : The highest harmonic order to retain in the Fourier series.
        %
        %   Output:
        %   - r ((2h+1)NM × 1 vector) : The Fourier representation of vect(A).
        %
        %   Example:
        %   % Compute Fourier-vectorized representation of A(t)
        %   A = PhasorArray(rand(4,4,11));
        %   FvTB_A = FvTB(A, 5);
        %
        %   % Compute using F_tb on the vectorized form
        %   FvTB_A_alt = F_tb(A{:}, 5); % Equivalent computation
        %
        %   See also: F_tb, vect, trunc, pad.
        if size(pA1,3)>2*h+1
            pA1=trunc(pA1,h);
        else
            pA1=pad(pA1,h-pA1.h);
        end
        d1=reshape(pA1,[],1,size(pA1,3));

        r=pvalue(permute(d1,[3 1 2]));
        r=r(:);
    end

    function [HpJ,JHm,Hp,Hm] = TBHankel(pA1,m)
        %TBHANKEL Compute Toeplitz Block Hankel matrices of order m.
        %
        %   This function computes the **Toeplitz Block Hankel matrices** (H) and the
        %   associated **J-Hankel matrices** (H_J) of the given matrix `A(t)` represented
        %   by the PhasorArray `pA1`.
        %
        %   **Definition:**
        %   - A **Hankel matrix** is symmetric with entries reflected across its antidiagonals.
        %   - A **J-Hankel matrix** is defined by pre- and post-multiplying the Hankel matrix
        %     with the flipping operator `J_m`, where `J_m` flips the diagonals of a matrix.
        %
        %   **Key Property:**
        %   - For a periodic matrix function A(t), the Toeplitz Block Hankel decomposition
        %     complements the Toeplitz Block transform (TB).
        %
        %   **Truncation:**
        %   - The matrix is truncated to order `2m` harmonics before performing the computation.
        %
        %   **Output Details:**
        %   - `HpJ`: Positive J-Hankel matrix for `A(t)`.
        %   - `JHm`: Negative J-Hankel matrix for `A(t)`.
        %   - `Hp`: Positive Hankel block matrix for `A(t)`.
        %   - `Hm`: Negative Hankel block matrix for `A(t)`.
        %
        %   **Dimensions:**
        %   - If `A` is an `N×N` matrix, then all outputs have dimensions `(2m+1)N × (2m+1)N`.
        %
        %   Syntax:
        %   [HpJ, JHm, Hp, Hm] = TBHANKEL(A, m)
        %       Computes the Toeplitz Block Hankel and J-Hankel matrices of order `m`.
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The input PhasorArray representing A(t).
        %   - m (integer) : The truncation order for the harmonics.
        %
        %   Output:
        %   - HpJ ((2m+1)N × (2m+1)N matrix) : Positive J-Hankel matrix.
        %   - JHm ((2m+1)N × (2m+1)N matrix) : Negative J-Hankel matrix.
        %   - Hp ((2m+1)N × (2m+1)N matrix) : Positive Hankel block matrix.
        %   - Hm ((2m+1)N × (2m+1)N matrix) : Negative Hankel block matrix.
        %
        %   Example:
        %   % Compute Hankel matrices for a given A(t)
        %   A = PhasorArray(rand(4,4,11));
        %   [HpJ, JHm, Hp, Hm] = TBHankel(A, 5);
        %
        %   See also: spTBHANKEL, BTHANKEL.
        [HpJ,JHm,Hp,Hm] = Array2TBHankel(pA1,2*m);
    end
    function [HpJ,JHm,Hp,Hm] = spTBHankel(pA1,m)
        %SPTBHANKEL Compute sparse Toeplitz Block Hankel matrices of order m.
        %
        %   This function computes the **sparse Toeplitz Block Hankel matrices** (H) and
        %   their associated **J-Hankel matrices** (H_J) for the given PhasorArray `pA1`.
        %
        %   **Sparse Computation:**
        %   - The matrices are stored and computed in **sparse format** for efficiency,
        %     especially for high-dimensional problems or large truncation orders.
        %
        %   **Output Details:** (same as `TBHANKEL`)
        %   - `HpJ`: Positive sparse J-Hankel matrix.
        %   - `JHm`: Negative sparse J-Hankel matrix.
        %   - `Hp`: Positive sparse Hankel block matrix.
        %   - `Hm`: Negative sparse Hankel block matrix.
        %
        %   **Syntax:**
        %   [HpJ, JHm, Hp, Hm] = SPTBHANKEL(A, m)
        %       Computes the sparse Toeplitz Block Hankel and J-Hankel matrices of order `m`.
        %
        %   See also: TBHANKEL, BTHANKEL.
        [HpJ,JHm,Hp,Hm] = spArray2TBHankel(pA1,2*m);
    end
    function [HpJ,JHm,Hp,Hm] = BTHankel(pA1,m)
        %BTHANKEL Compute Block Toeplitz Hankel matrices of order m.
        %
        %   This function computes the **Block Toeplitz Hankel matrices** (H) and the
        %   associated **J-Hankel matrices** (H_J) for the given PhasorArray `pA1`.
        %
        %   **Difference with TBHankel:**
        %   - While `TBHankel` corresponds to **Toeplitz Block (TB)** structures, `BTHankel`
        %     corresponds to **Block Toeplitz (BT)** structures.
        %
        %   **Output Details:** (same as `TBHANKEL`)
        %   - `HpJ`: Positive J-Hankel matrix.
        %   - `JHm`: Negative J-Hankel matrix.
        %   - `Hp`: Positive Hankel block matrix.
        %   - `Hm`: Negative Hankel block matrix.
        %
        %   **Syntax:**
        %   [HpJ, JHm, Hp, Hm] = BTHANKEL(A, m)
        %       Computes the Block Toeplitz Hankel and J-Hankel matrices of order `m`.
        %
        %   See also: SPBTHANKEL, TBHANKEL.
        [HpJ,JHm,Hp,Hm] = Array2BTHankel(pA1,2*m);
    end
    function [HpJ,JHm,Hp,Hm] = spBTHankel(pA1,m)
        %SPBTHANKEL Compute sparse Block Toeplitz Hankel matrices of order m.
        %
        %   This function computes the **sparse Block Toeplitz Hankel matrices** (H) and
        %   the associated **J-Hankel matrices** (H_J) for the given PhasorArray `pA1`.
        %
        %   **Sparse Computation:** (same as `spTBHANKEL`)
        %   - The matrices are stored and computed in **sparse format** for efficiency.
        %
        %   **Syntax:**
        %   [HpJ, JHm, Hp, Hm] = SPBTHANKEL(A, m)
        %       Computes the sparse Block Toeplitz Hankel and J-Hankel matrices of order `m`.
        %
        %   See also: BTHANKEL, TBHANKEL.
        [HpJ,JHm,Hp,Hm] = spArray2BTHankel(pA1,2*m);
    end

    function AB_TB = TBmtimes(pA1,pA2,h)
        % TBMTIMES Compute the Toeplitz Block (TB) matrix of the product A(t)B(t).
        %
        %   TBMTIMES(pA1, pA2, h) computes the Toeplitz Block matrix of order `h`
        %   that represents the product of two periodic matrix functions A(t) and B(t).
        %
        %   Inputs:
        %     pA1  - (PhasorArray) The first input PhasorArray representing A(t).
        %     pA2  - (PhasorArray) The second input PhasorArray representing B(t).
        %     h   - (integer) The truncation order for the harmonics.
        %
        %   Outputs:
        %     AB_TB - ((2h+1)N × (2h+1)N matrix) The Toeplitz Block matrix for A(t)B(t).
        %
        %   Behavior:
        %     - Uses the Hankel matrices and TB representation of A(t) and B(t).
        %     - The resulting TB matrix satisfies:
        %         T_tb(A * B) = T_tb(A) * T_tb(B) + H_p(A) * H_m(B) + JH_m(A) * JH_p(B).
        %
        %   Output Dimensions:
        %     - If A and B are `N×N` matrices and the truncation order is `h`, then:
        %         - `AB_TB` is a `(2h+1)N × (2h+1)N` matrix.
        %
        %   Example Usage:
        %     % Compute TB matrix of the product A(t)B(t)
        %     A = PhasorArray(rand(4,4,11));
        %     B = PhasorArray(rand(4,4,11));
        %     AB_TB = TBmtimes(A, B, 5);
        %
        %   See also: spTBmtimes, TB, TBHankel.
        [~,JHm,Hp,~] = TBHankel(pA1,h);
        [HpJ2,~,~,Hm2] = TBHankel(pA2,h);
        o1TB=pA1.T_tb(h);
        o2TB=pA2.T_tb(h);
        AB_TB=o1TB*o2TB+Hp*Hm2+JHm*HpJ2;
    end
    function AB_TB = spTBmtimes(pA1,pA2,h)
        % SPTBMTIMES Compute the sparse Toeplitz Block (TB) matrix of the product A(t)B(t).
        %
        %   SPTBMTIMES(pA1, pA2, h) computes the sparse Toeplitz Block matrix of order `h`
        %   that represents the product of two periodic matrix functions A(t) and B(t).
        %
        %   Inputs:
        %     pA1  - (PhasorArray) The first input PhasorArray representing A(t).
        %     pA2  - (PhasorArray) The second input PhasorArray representing B(t).
        %     h   - (integer) The truncation order for the harmonics.
        %
        %   Outputs:
        %     AB_TB - ((2h+1)N × (2h+1)N sparse matrix) Sparse TB matrix for A(t)B(t).
        %
        %   Behavior:
        %     - Uses sparse versions of Hankel matrices and TB representation of A(t) and B(t).
        %     - The result satisfies:
        %         T_tb(A * B) = T_tb(A) * T_tb(B) + H_p(A) * H_m(B) + JH_m(A) * JH_p(B).
        %
        %   Output Dimensions:
        %     - If A and B are `N×N` matrices and the truncation order is `h`, then:
        %         - `AB_TB` is a `(2h+1)N × (2h+1)N` sparse matrix.
        %
        %   Example Usage:
        %     % Compute sparse TB matrix of the product A(t)B(t)
        %     A = PhasorArray(rand(4,4,11));
        %     B = PhasorArray(rand(4,4,11));
        %     AB_TB = spTBmtimes(A, B, 5);
        %
        %   See also: TBmtimes, spTB, spTBHankel.
        [~,JHm,Hp,~] = spTBHankel(pA1,h);
        [HpJ2,~,~,Hm2] = spTBHankel(pA2,h);
        o1TB=pA1.spTB(h);
        o2TB=pA2.spTB(h);
        AB_TB=o1TB*o2TB+Hp*Hm2+JHm*HpJ2;
    end

    function r= F_tb(pA1,m)
        %F_tb Compute the Fourier representation of A in a form compatible with T_tb(A, m).
        %
        %   This function computes the **Fourier representation** of A(t) up to order `m`,
        %   but instead of forming a **square Toeplitz Block matrix** (T_tb(A, m)), it **stacks
        %   the phasors of each element of A into a structured vectorized form**.
        %
        %   **Definition:**
        %   - Let `a(t)` be a scalar function with Fourier coefficients `(a_k)_{k∈ℤ}`.
        %   - Then, `F_tb(a, m) = [a_(-m); a_(-m+1); ... a_0; a_1; ... a_m]`.
        %   - If `A(t)` is an `N×M` matrix, its Fourier representation is given by:
        %       F_tb(A, m) = [F_tb(A_11), F_tb(A_12); F_tb(A_21), F_tb(A_22)].
        %
        %   **Key Property (TB Compatibility):**
        %   If `y(t) = A(t)x(t)`, then:
        %       F_TB(y) = T_tb(A, m) ⋅ F_TB(x),
        %   where `T_tb(A, m)` is the **Toeplitz Block representation** of `A` (see `TB` method).
        %
        %   **Dimension of the Output:**
        %   - If `A` is an `N×M` matrix, then `F_tb(A, m)` is a `((2m+1)N) × M` matrix.
        %
        %   Syntax:
        %   r = F_tb(A, m)
        %       Computes the Fourier representation of `A` up to order `m`, compatible with T_tb(A, m).
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The input PhasorArray representing A(t).
        %   - m (integer, optional) : The highest harmonic order to retain. Default: `A.h`.
        %
        %   Output:
        %   - r ((2m+1)N × M matrix) : The Fourier representation of A in a form compatible with T_tb(A, m).
        %
        %   Example:
        %   % Compute Fourier representation of A in TB form
        %   A = PhasorArray(rand(4,4,11));
        %   F_tb_A = F_tb(A, 5);
        %
        %   See also: F_bt, TB.
        arguments
            pA1
            m=(size(pA1,3)-1)/2;
        end
        nho=(size(pA1,3)-1)/2;
        if nho<m
            toto=value(pA1.pad(m-nho));
        elseif nho>m
            toto=pvalue(trunc(pA1,m));
        else
            toto=pvalue(pA1);
        end
        titi=permute(toto,[3 1 2]);
        r=reshape(titi,[],size(pA1,2),1);
    end
    function r= F_bt(pA1,m)
        %F_bt Compute the Fourier representation of A in a form compatible with T_bt(A, m).
        %
        %   This function computes the **Fourier representation** of A(t) up to order `m`,
        %   but instead of forming a **square Block Toeplitz matrix** (T_bt(A, m)), it **stacks
        %   the phasors of each element of A into a structured vectorized form**.
        %
        %   **Definition:**
        %   - Let `a(t)` be a scalar function with Fourier coefficients `(a_k)_{k∈ℤ}`.
        %   - Then, `F_bt(a, m) = [a_(-m); a_(-m+1); ... a_0; a_1; ... a_m]`.
        %   - If `A(t)` is an `N×M` matrix, its Fourier representation is given by:
        %       F_bt(A, m) = [F_bt(A_11), F_bt(A_12); F_bt(A_21), F_bt(A_22)].
        %
        %   **Key Property (BT Compatibility):**
        %   If `y(t) = A(t)x(t)`, then:
        %       F_BT(y) = T_bt(A, m) ⋅ F_BT(x),
        %   where `T_bt(A, m)` is the **Block Toeplitz representation** of `A` (see `BT` method).
        %
        %   **Dimension of the Output:**
        %   - If `A` is an `N×M` matrix, then `F_bt(A, m)` is an `N × (2m+1)M` matrix.
        %
        %   Syntax:
        %   r = F_bt(A, m)
        %       Computes the Fourier representation of `A` up to order `m`, compatible with T_bt(A, m).
        %
        %   Input Arguments:
        %   - pA1 (PhasorArray) : The input PhasorArray representing A(t).
        %   - m (integer, optional) : The highest harmonic order to retain. Default: `A.h`.
        %
        %   Output:
        %   - r (N × (2m+1)M matrix) : The Fourier representation of A in a form compatible with T_bt(A, m).
        %
        %   Example:
        %   % Compute Fourier representation of A in BT form
        %   A = PhasorArray(rand(4,4,11));
        %   F_bt_A = F_bt(A, 5);
        %
        %   See also: F_tb, BT.
        arguments
            pA1
            m=(size(pA1,3)-1)/2;
        end
        nho=(size(pA1,3)-1)/2;
        if nho<m
            toto=pvalue(pA1.pad(m-nho));
        elseif nho>m
            toto=pvalue(trunc(pA1,m));
        else
            toto=pvalue(pA1);
        end
        titi=permute(toto,[1 3 2]);
        r=reshape(titi,[],size(pA1,2),1);
    end

    function r = spFTB(pA1, m)
        %SPFTB Sparse Fourier-Toeplitz-Block representation.
        %
        %   r = spFTB(pA1, m) returns a sparse Fourier representation
        %   compatible with spTB(pA1, m).
        %
        %   See also: F_tb, spTB.
        arguments
            pA1
            m = (size(pA1,3)-1)/2;
        end
        r = pA1.F_tb(m);
        if ~isspecial(pA1)
            r = sparse(r);
        end
    end

    function r = spFBT(pA1, m)
        %SPFBT Sparse Fourier-Block-Toeplitz representation.
        %
        %   r = spFBT(pA1, m) returns a sparse Fourier representation
        %   compatible with spBT(pA1, m).
        %
        %   See also: F_bt, spBT.
        arguments
            pA1
            m = (size(pA1,3)-1)/2;
        end
        r = pA1.F_bt(m);
        if ~isspecial(pA1)
            r = sparse(r);
        end
    end

    function r= HmqNEig(pA1,h,T,bandlimit)
        % HMQNEIG Compute the eigenvalues of T_tb(pA1, h) - N_tb(pA1, h, T).
        %
        %   HMQNEIG(pA1, h, T, bandlimit) computes the eigenvalues of the difference between
        %   the Toeplitz Block matrix `T_tb(pA1, h)` and the non-periodic Toeplitz approximation
        %   `N_tb(pA1, h, T)`, capturing spectral properties of `pA1`.
        %
        %   Inputs:
        %     pA1        - (PhasorArray) The input PhasorArray representing a periodic matrix function.
        %     h         - (integer, optional) The truncation order for the harmonics.
        %                   - Default: `2*pA1.h` (twice the highest harmonic stored in `pA1`).
        %     T         - (double) The period of the PhasorArray.
        %     bandlimit - (char) Specifies band-limiting of eigenvalues:
        %                   - 'none' (default): No band-limiting.
        %                   - 'fundamental': Retains eigenvalues where |imag(r)| < π / |T|.
        %
        %   Outputs:
        %     r - (vector) The eigenvalues of the matrix `T_tb(pA1, h) - N_tb(pA1, h, T)`.
        %
        %   Behavior:
        %     - Computes eigenvalues of the difference matrix `Hmq = T_tb(pA1, h) - N_tb(pA1, h, T)`.
        %     - Band-limiting removes eigenvalues outside the specified range.
        %
        %   Example Usage:
        %     % Compute eigenvalues with fundamental band-limiting
        %     A = PhasorArray(rand(4,4,11));
        %     r = HmqNEig(A, 10, 2, 'fundamental');
        %
        %     % Compute eigenvalues without band-limiting
        %     r = HmqNEig(A, [], 1, 'none');
        %
        %   See also: TB, N_tb, eig.
        arguments
            pA1
            h=20
            T=2*pi
            bandlimit {mustBeMember(bandlimit,{'none','fundamental'})}='none'
        end
        if isempty(h)
            h=2*pA1.h;
        end
        r1=pA1.T_tb(h)-N_tb(pA1.value,h,T);
        r=eig(r1);
        switch bandlimit
            case 'fundamental'
                r=r(abs(imag(r))<pi/abs(T));
            case 'none'
            otherwise
                r=r(abs(imag(r))<=bandlimit);

        end
    end


    function r= HmqEig(pA1,h)
        % HMQEIG Compute the eigenvalues of the Toeplitz Block (TB) matrix of A(t).
        %
        %   HMQEIG(pA1, h) computes the eigenvalues of the Toeplitz Block (TB) matrix
        %   representation of the periodic matrix function `A(t)`, capturing its spectral properties.
        %
        %   Inputs:
        %     pA1  - (PhasorArray) The input PhasorArray representing a periodic matrix function.
        %     h   - (integer, optional) The truncation order for the harmonics.
        %            - Default: `pA1.h` (the highest harmonic stored in `pA1`).
        %
        %   Outputs:
        %     r   - (vector) The eigenvalues of the Toeplitz Block matrix `T_tb(pA1, h)`.
        %
        %   Behavior:
        %     - The matrix is truncated to order `h` to form `T_tb(pA1, h)`.
        %     - If `h` is not specified, the function uses the highest harmonic available in `pA1`.
        %
        %   Output Dimensions:
        %     - If `A(t)` is an `N×N` matrix and the truncation order is `h`, then:
        %         - `T_tb(pA1, h)` is a `(2h+1)N × (2h+1)N` matrix.
        %         - `r` is a vector of eigenvalues of length `(2h+1)N`.
        %
        %   Example Usage:
        %     % Compute eigenvalues of T_tb(A, h) with a specific truncation order
        %     A = PhasorArray(rand(4,4,11));
        %     r = HmqEig(A, 10);
        %
        %     % Compute eigenvalues using the default harmonic order
        %     r = HmqEig(A);
        %
        %   See also: TB, eig.
        arguments
            pA1
            h=pA1.h
        end
        r1=pA1.T_tb(h);
        r=eig(r1);
    end

    function r=evalp(pA1,angle,nvp)
        % EVALP Evaluate a periodic matrix at a given phase angle.
        %
        %   EVALP(pA1, angle, nvp) computes the value of a periodic matrix `A` at a
        %   specified phase angle `θ` instead of time. This function is useful for
        %   evaluating periodic signals in the frequency domain.
        %
        %   Syntax:
        %     r = EVALP(A, angle)
        %     r = EVALP(A, angle, Name, Value)
        %
        %   Inputs:
        %     pA1    - (PhasorArray) The periodic matrix represented as a **phasor array**.
        %     angle - (double, vector) The phase angle(s) at which to evaluate `A(θ)`.
        %
        %   Name-Value Pair Arguments:
        %     'forceReal'    - (logical) If `true`, forces real-valued output. Default: `false`.
        %     'checkReal'    - (logical) If `true`, checks whether the result is real-valued. Default: `false`.
        %     'checkRealTol' - (double) Tolerance for checking real-valued results. Default: `1e-8`.
        %
        %   Outputs:
        %     r - (matrix) Evaluated periodic matrix at the given **phase angle** `θ`.
        %
        %   Behavior:
        %     - This function evaluates `A(θ)`, where `θ` represents the **phase angle** of the periodic function.
        %     - The function internally calls `PhasorArray2time` with `T = 2π`, assuming **one full cycle** of the periodic function.
        %     - If `pA1` is detected as real-valued, `forceReal` is automatically set to `true`.
        %
        %   Example Usage:
        %     % Evaluate a phasor array at a phase angle of π/4
        %     result = evalp(A, pi/4);
        %
        %     % Force real-valued output
        %     result = evalp(A, pi/4, 'forceReal', true);
        %
        %   See also: PHASORARRAY2TIME, EVALTIME.
        arguments
            pA1
            angle
            nvp.forceReal logical    = false
            nvp.checkReal logical    = false
            nvp.checkRealTol         = 1e-8
            nvp.squeeze              = false
        end

        if isreal(pA1)
            nvp.forceReal = true;
        end

        nvp = namedargs2cell(nvp);
        %evaluate periodic matrix A for an angle argument instead of
        %time
        r=PhasorArray2time(pA1,2*pi,angle,nvp{:});
    end


    function r = mreal(pA1)
        %MREAL Compute the real part of a PhasorArray in the time domain.
        %   r = MREAL(pA1) returns the real part of the PhasorArray pA1 in the time domain.
        %   The real part R(t) is such that A(t) = R(t) + i*I(t), where I(t) is the imaginary part.
        %   Both R(t) and I(t) are real-valued.
        %
        %   INPUT:
        %       pA1 - The PhasorArray object.
        %
        %   OUTPUT:
        %       r  - The real part of the PhasorArray in the time domain.
        %
        %   Example:
        %       A = PhasorArray.random(3, 3, 5);
        %       R = mreal(A);
        %
        %   See also: mimag, mconj
        dval = pvalue(pA1);
        % real/imag pulled out before combining, and * instead of /: YALMIP
        % supports neither on a derived ndsdpvar. Bit-identical by linearity.
        dr = real(dval);   di = imag(dval);
        r  = (dr + flip(dr,3))*(1/2) + 1i * (di - flip(di,3))*(1/2);
        r = PhasorArray(r);
    end

    function r = mimag(pA1)
        %MIMAG Compute the imaginary part of a PhasorArray in the time domain.
        %   r = MIMAG(pA1) returns the imaginary part of the PhasorArray pA1 in the time domain.
        %   The imaginary part I(t) is such that A(t) = R(t) + i*I(t), where R(t) is the real part.
        %   Both R(t) and I(t) are real-valued.
        %
        %   INPUT:
        %       pA1 - The PhasorArray object.
        %
        %   OUTPUT:
        %       r  - The imaginary part of the PhasorArray in the time domain.
        %
        %   Example:
        %       A = PhasorArray.random(3, 3, 5);
        %       I = mimag(A);
        %
        %   See also: mreal, mconj
        dval = pvalue(pA1);
        % Same rewrites as mreal (1/1i == -1i exactly, so *(-1i) replaces /1i).
        dr = real(dval);   di = imag(dval);
        r  = (dr - flip(dr,3))*(1/2)*(-1i) + 1i * (di + flip(di,3))*(1/2)*(-1i);
        r = PhasorArray(r);
    end

    function r = mconj(pA1)
        %MCONJ Compute the complex conjugate of a PhasorArray in the time domain.
        %   r = MCONJ(pA1) returns the complex conjugate of the PhasorArray pA1 in the time domain.
        %   The complex conjugate is computed as R(t) - i*I(t), where R(t) is the real part and I(t) is the imaginary part.
        %
        %   INPUT:
        %       pA1 - The PhasorArray object.
        %
        %   OUTPUT:
        %       r  - The complex conjugate of the PhasorArray in the time domain.
        %
        %   Example:
        %       A = PhasorArray.random(3, 3, 5);
        %       C = mconj(A);
        %
        %   See also: mreal, mimag
        r = mreal(pA1)-1i*mimag(pA1);
    end

    function pA1=confirm_reality(pA1)
        %  depreciated
        %
        % see also MREAL.

        pos_part=pA1{:,:,1:pA1.h};
        neg_part=pA1{:,:,-1:-1:(-pA1.h)};
        z_part=pA1.phas(0);

        pA1{:,:,1:pA1.h}=(pos_part+conj(neg_part))/2;
        pA1{:,:,0}=real(z_part);
        pA1{:,:,-1:-1:(-pA1.h)}=conj(pos_part+conj(neg_part))/2;
    end

    function r = pageconj(pA1)
        %phasor conjugate
        r=(conj(pvalue(pA1)));
    end
    function r = pagereal(pA1)
        %phasor real part
        r=(real(pvalue(pA1)));
    end
    function r = pageimag(pA1)
        %phasor imag part
        r=(imag(pvalue(pA1)));
    end
    function r = pageabs(pA1)
        %phasor absolute value
        r=(abs(pvalue(pA1)));
    end


    function r = isvector(pA1)
        % ISVECTOR True if A(t) is a vector (row or column).
        r=isvector(pA1.phas(0));
    end
    function r = isscalar(pA1)
        % ISSCALAR True if A(t) is a scalar.
        r=isscalar(pA1.phas(0));
    end
    function r = isrow(pA1)
        % ISROW True if A(t) is a row vector.
        r=isrow(pA1.phas(0));
    end
    function r = iscolumn(pA1)
        % ISCOLUMN True if A(t) is a column vector.
        r=iscolumn(pA1.phas(0));
    end
    function r = ismatrix(pA1)
        % ISMATRIX True if A(t) is a matrix.
        r=ismatrix(pA1.phas(0));
    end
    function r = isnumeric(pA1)
        % ISNUMERIC True if A(t) contains numeric values.
        r=isnumeric(pA1.value);
    end
    function r = islogical(pA1)
        % ISLOGICAL True if A(t) contains logical values.
        r=islogical(pA1.value);
    end
    function r = isempty(pA1)
        % ISEMPTY True if A(t) is empty.
        r=isempty(pA1.value);
    end
    function r = issquare(pA1)
        % ISSQUARE True if A(t) is a square matrix.
        r = (size(pA1,1)==size(pA1,2));
    end

    function [r,R] = issymmetric(pA1,nvp)
        % ISSYMMETRIC Check if the PhasorArray object is symmetric.
        %
        %   [r, R] = ISSYMMETRIC(pA1, nvp) checks if the PhasorArray object pA1
        %   is symmetric based on the specified arguments.
        %
        %   Input arguments:
        %       pA1 - The PhasorArray object to be checked.
        %       nvp - Name-Value arguments:
        %           skewOption - A string specifying the type of symmetry to check.
        %                         It can be either 'nonskew' or 'skew'. Default is 'nonskew'.
        %           tol - A tolerance value for the symmetry check. Default is 0.
        %
        %   Output arguments:
        %       r - A logical value indicating if the entire PhasorArray object is symmetric.
        %       R - A logical array indicating the symmetry of each slice of the PhasorArray object.
        %
        %   Example:
        %       [r, R] = issymmetric(pA1, 'skewOption', 'nonskew', 'tol', 1e-6);
        %
        %   See also: arrayfun, norm, nnz
        arguments
            pA1
            nvp.skewOption {mustBeMember(nvp.skewOption,{'nonskew','skew'})} = 'nonskew'
            nvp.tol = 0
        end
        if nvp.tol == 0
            R = arrayfun(@(ii) issymmetric(pA1(:,:,ii),nvp.skewOption),1:(2*pA1.h+1));
        else
            tol=abs(nvp.tol);
            switch nvp.skewOption
                case 'nonskew'
                    R = arrayfun(@(ii) norm(pA1(:,:,ii)-pA1(:,:,ii).',"inf")/norm(pA1(:,:,ii),"inf")<tol,1:(2*pA1.h+1));
                otherwise
                    R = arrayfun(@(ii) norm(pA1(:,:,ii)+pA1(:,:,ii).',"inf")/norm(pA1(:,:,ii),"inf")<tol,1:(2*pA1.h+1));
            end
        end
        r = ((nnz(R))==(2*pA1.h+1));
    end

    function [r,R] = ishermitian(pA1,nvp)
        % ISHERMITIAN Check if the PhasorArray object is Hermitian.
        %
        %   [r, R] = ISHERMITIAN(pA1, nvp) checks if the PhasorArray object pA1
        %   is Hermitian. The function returns a logical scalar r indicating
        %   if all slices of the PhasorArray are Hermitian, and a logical array R
        %   indicating if each individual slice is Hermitian.
        %
        %   Input arguments:
        %       pA1 - PhasorArray object to be checked.
        %       nvp - Name-Value arguments:
        %           skewOption - (optional) Specifies whether to check for
        %                         'nonskew' (default) or 'skew' Hermitian.
        %                         Must be one of {'nonskew', 'skew'}.
        %           tol - (optional) Tolerance for numerical comparison. Default is 0.
        %
        %   Output arguments:
        %       r - Logical scalar indicating if all slices of the PhasorArray
        %           are Hermitian.
        %       R - Logical array indicating if each individual slice of the
        %           PhasorArray is Hermitian.
        arguments
            pA1
            nvp.skewOption {mustBeMember(nvp.skewOption,{'nonskew','skew'})} = 'nonskew'
            nvp.tol = 0
        end
        if nvp.tol == 0
            R = arrayfun(@(ii) ishermitian(pA1(:,:,ii),nvp.skewOption),1:(2*pA1.h+1));
        else
            tol=abs(nvp.tol);
            switch nvp.skewOption
                case 'nonskew'
                    R = arrayfun(@(ii) norm(pA1(:,:,ii)-pA1(:,:,ii)',"inf")/norm(pA1(:,:,ii),"inf")<tol,1:(2*pA1.h+1));
                otherwise
                    R = arrayfun(@(ii) norm(pA1(:,:,ii)+pA1(:,:,ii)',"inf")/norm(pA1(:,:,ii),"inf")<tol,1:(2*pA1.h+1));
            end
        end
        r = ((nnz(R))==(2*pA1.h+1));
    end

    function [r,R,tolmin] = isreal(pA1,tol)
        % ISREAL Check if the imaginary part of the time realization of
        % PhasorArray is negligible compared to a given tolerance.
        %
        % Syntax:
        %   [r, R] = isreal(pA1)
        %   [r, R] = isreal(pA1, tol)
        %
        % Description:
        %   This function checks if the imaginary part of the time realization
        %   of a PhasorArray object is negligible compared to a specified tolerance.
        %   It returns a boolean value 'r' indicating whether the imaginary part
        %   is negligible, and a logical array 'R' indicating the same for each
        %   element of the PhasorArray.
        %
        % Inputs:
        %   pA1  - PhasorArray object to be checked.
        %   tol - (Optional) Tolerance value for comparison. If not provided,
        %         the function uses the default machine epsilon.
        %
        % Outputs:
        %   r   - Boolean value (true/false) indicating if the imaginary part
        %         is negligible for the entire PhasorArray.
        %   R   - Logical array indicating if the imaginary part is negligible
        %         for each element of the PhasorArray.
        %
        % Example:
        %   [r, R] = isreal(pA1)
        %   [r, R] = isreal(pA1, 1e-6)
        %
        % See also:
        %   mreal, ismembertol, abs
        arguments
            pA1
            tol=[]
        end
        [r, R] = projectionMatch(pA1, @mreal, +1, tol, 'isreal');
        if nargout > 2
            if r == 0 && isspecial(pA1)
                tolmin = NaN;   % undetermined: bisecting on it would not converge
            else
                tolmin = tolReal(pA1);
            end
        end
    end
    function tol = tolZero(pA1, tolstart, tolTol)
        % TOLZERO Determines the lowest tolerance for which the input is considered zero.
        %
        % Syntax: tol = tolZero(pA1, tolstart, tolTol)
        %
        % Inputs:
        %    pA1 - The input PhasorArray object to be checked for zero values.
        %    tolstart - (Optional) The starting tolerance value for the binary search. Default is 100.
        %    tolTol - (Optional) The tolerance for the binary search convergence. Default is 1e-20.
        %
        % Outputs:
        %    tol - The lowest tolerance for which the input is considered zero through function iszero(pA1, tol).
        %
        % Example:
        %    tol = tolZero(pA1, 100, 1e-20);
        %
        % See also: iszero
        arguments
            pA1
            tolstart = 100
            tolTol = 1e-6
        end
        % Initialize variables
        low = 0;
        high = tolstart;

        % Binary search for minimum tolerance
        while (high - low) > tolTol
            mid = (low + high) / 2;
            if iszero(pA1, mid)
                high = mid;
            else
                low = mid;
            end
        end

        % Return the minimum tolerance
        tol = high;
    end


    function tol = tolReal(pA1,tolstart,tolTol)
        % tolReal - Determines the lowest tolerance for which the input is real.
        %
        % Syntax: tol = tolReal(pA1, tolstart, tolTol)
        %
        % Inputs:
        %    pA1 - The input object to be checked for reality.
        %    tolstart - (Optional) The starting tolerance value for the binary search. Default is 100.
        %    tolTol - (Optional) The tolerance for the binary search convergence. Default is 1e-20.
        %
        % Outputs:
        %    tol - The lowest tolerance for which the input is real through function isreal(pA1, tol).
        %
        % Example:
        %    tol = tolReal(pA1, 100, 1e-20);
        %
        % See also: isreal
        arguments
            pA1
            tolstart=100
            tolTol=1e-20;
        end
        % Initialize variables
        low = 0;
        high = tolstart;

        % Binary search for minimum tolerance
        while (high - low) > tolTol
            mid = (low + high) / 2;
            if isreal(pA1, mid)
                high = mid;
            else
                low = mid;
            end
        end

        % Return the minimum tolerance
        tol = high;
    end

    function [r,R] = isimag(pA1,tol)
        % ISIMAG Check if the real part of the time realization of a PhasorArray is negligible.
        %
        %   [r, R] = ISIMAG(pA1, tol) checks whether the real part of the time-domain
        %   realization of the PhasorArray `pA1` is negligible relative to a given tolerance.
        %
        %   Inputs:
        %       pA1  - (PhasorArray) The PhasorArray object to be checked.
        %       tol - (double, optional) Tolerance threshold for comparison.
        %             - Default: `[]`, meaning it uses machine precision.
        %
        %   Outputs:
        %       r - (logical) True if the entire PhasorArray is purely imaginary.
        %       R - (logical array) Element-wise result indicating if each entry is purely imaginary.
        %
        %   Example:
        %       [r, R] = isimag(pA1, 1e-6);
        %
        %   See also: mimag, isreal, abs
        arguments
            pA1
            tol=[]
        end
        [r, R] = projectionMatch(pA1, @mimag, -1, tol, 'isimag');
    end
    function [r,R] = iscomplex(pA1,tol)
        % ISCOMPLEX Check if the PhasorArray contains significant imaginary components.
        %
        %   [r, R] = ISCOMPLEX(pA1, tol) determines whether the PhasorArray `pA1`
        %   contains significant imaginary components relative to a given tolerance.
        %
        %   Inputs:
        %       pA1  - (PhasorArray) The PhasorArray object to be checked.
        %       tol - (double, optional) Tolerance threshold for comparison.
        %             - Default: `eps` (machine precision).
        %
        %   Outputs:
        %       r - (logical) True if the entire PhasorArray has significant imaginary components.
        %       R - (logical array) Element-wise result indicating if each entry is complex.
        %
        %   Example:
        %       [r, R] = iscomplex(pA1, 1e-6);
        %
        %   See also: isreal, isimag, abs
        arguments
            pA1
            tol=eps
        end
        [r,R]=isreal(pA1,tol);
        r=~r;
        R=~R;
    end
    function [r,R] = iszero(pA1,tol)
        % ISZERO Check if the PhasorArray is numerically zero.
        %
        %   [r, R] = ISZERO(pA1, tol) determines whether the PhasorArray `pA1`
        %   is numerically zero, based on an absolute tolerance.
        %
        %   Inputs:
        %       pA1  - (PhasorArray) The PhasorArray object to be checked.
        %       tol - (double, optional) Tolerance threshold for comparison.
        %             - Default: `0`, meaning exact zero comparison.
        %
        %   Outputs:
        %       r - (logical) True if all entries in the PhasorArray are within the tolerance of zero.
        %       R - (logical array) Element-wise result indicating if each entry is within tolerance.
        %
        %   Example:
        %       [r, R] = iszero(pA1, 1e-6);
        %
        %   See also: pagenorm, isreal, iscomplex
        arguments
            pA1
            tol = 0
        end
        if ~isspecial(pA1)
            R = pageInfNorm(pA1.value)<=tol;
        elseif isa(pA1.value, 'sdpvar') || isa(pA1.value, 'ndsdpvar')
            % Undecidable on a decision variable: the value does not exist until
            % the problem is solved (YALMIP turns `x == 0` into a constraint).
            error('PhasorArray:iszero:decisionVariable', ...
                ['iszero is undecidable on an sdpvar/ndsdpvar-backed PhasorArray. ', ...
                 'Apply it to the solved array, e.g. iszero(PhasorArray(value(x))).'])
        else
            % Symbolic: pagenorm needs floats, so compare slice-wise instead.
            R = false(1, 1, 2*pA1.h+1);
            for k = 1:numel(R)
                R(k) = symbolicSliceIsZero(pA1(:, :, k), tol);
            end
        end
        r = ((nnz(R))==(2*pA1.h+1));
    end

    function r = isspecial(pA1)
        %ISSPECIAL True when the coefficients are symbolic or decision variables.
        %
        %   True for sym, sdpvar and ndsdpvar backings. Numeric predicates such as
        %   iszero cannot be evaluated on those, so guard with this first.
        r = (isa(pA1.value,"sym") || isa(pA1.value,"ndsdpvar") || isa(pA1.value,"sdpvar"));
    end

    function r = ImagRealForm(pA1)
        % IMAGREALFORM Convert PhasorArray to Imaginary-Real form.
        %
        %   r = IMAGREALFORM(pA1) converts the PhasorArray `pA1` into the
        %   Imaginary-Real representation, ensuring:
        %       - **Negative harmonics (`-k`) store the imaginary part (left side).**
        %       - **Positive harmonics (`+k`) store the real part (right side).**
        %       - The order `k` and `-k` remain **symmetric** around the DC component.
        %
        %   This format assumes that `pA1` represents a **real-valued periodic matrix**,
        %   meaning `X_-k = conj(X_k)`.
        %
        %   Input:
        %       pA1 - (PhasorArray) The input PhasorArray object.
        %
        %   Output:
        %       r - (array) The Imaginary-Real form representation:
        %           - Dimensions: `(m × n × (2h+1))`
        %           - Order `-k` maps to `(:,:,h-k+1)` (imaginary part, negated).
        %           - Order `+k` maps to `(:,:,h+k+1)` (real part).
        %
        %   Example:
        %       r = ImagRealForm(pA1);
        %
        %   See also: RealImagForm, SinCosForm
        h=pA1.h;
        pA1=pA1.value;
        indices = repmat({':'}, 1, ndims(pA1)-3);
        r=cat(3,-imag(pA1(:,:,1:h,indices{:})),real(pA1(:,:,(h+1):end,indices{:})));
    end
    function r = RealImagForm(pA1)
        % REALIMAGFORM Convert PhasorArray to Real-Imaginary form.
        %
        %   r = REALIMAGFORM(pA1) converts the PhasorArray `pA1` into the
        %   Real-Imaginary representation, ensuring:
        %       - **Negative harmonics (`-k`) store the real part (left side).**
        %       - **Positive harmonics (`+k`) store the imaginary part (right side).**
        %       - The order `k` and `-k` remain **symmetric** around the DC component.
        %
        %   This format assumes that `pA1` represents a **real-valued periodic matrix**,
        %   meaning `X_-k = conj(X_k)`.
        %
        %   Input:
        %       pA1 - (PhasorArray) The input PhasorArray object.
        %
        %   Output:
        %       r - (array) The Real-Imaginary form representation:
        %           - Dimensions: `(m × n × (2h+1))`
        %           - Order `-k` maps to `(:,:,h-k+1)` (real part).
        %           - Order `+k` maps to `(:,:,h+k+1)` (imaginary part).
        %
        %   Example:
        %       r = RealImagForm(pA1);
        %
        %   See also: ImagRealForm, SinCosForm
        r = flip(ImagRealForm(pA1),3);
    end
    function r = SinCosForm(pA1,isReal,realTol)
        % SINCOSFORM Convert PhasorArray to Sine-Cosine form.
        %
        %   r = SINCOSFORM(pA1, isReal, realTol) converts the PhasorArray `pA1`
        %   into the **Sine-Cosine representation**, ensuring:
        %       - **Negative harmonics (`-k`) store sine coefficients (left side).**
        %       - **Positive harmonics (`+k`) store cosine coefficients (right side).**
        %       - The order `k` and `-k` remain **symmetric** around the DC component.
        %
        %   Inputs:
        %       pA1      - (PhasorArray) The input PhasorArray object.
        %       isReal  - (logical, optional) If `true`, forces real-valued output.
        %                 Default: `isreal(pA1)`.
        %       realTol - (double, optional) Tolerance for enforcing real values.
        %                 Default: `1e-14`.
        %
        %   Output:
        %       r - (array) The Sine-Cosine form representation:
        %           - Dimensions: `(m × n × (2h+1))`
        %           - Order `-k` maps to `(:,:,h-k+1)` (sine coefficient, negated).
        %           - Order `+k` maps to `(:,:,h+k+1)` (cosine coefficient).
        %
        %   Example:
        %       r = SinCosForm(pA1, true, 1e-14);
        %
        %   See also: CosSinForm, ImagRealForm
        arguments
            pA1
            isReal logical = isreal(pA1)
            realTol = 1e-14
        end
        %normal sin / cos form
        h=pA1.h;
        pA1=pA1.value;
        indices = repmat({':'}, 1, ndims(pA1)-3);
        r=cat(3,1i*(flip(pA1(:,:,h+2:end,indices{:}),3)-pA1(:,:,1:h,indices{:})),pA1(:,:,h+1,indices{:}),(pA1(:,:,h+2:end,indices{:})+flip(pA1(:,:,1:h,indices{:}),3)));
        if isReal
            er = real(r-real(r));
            if norm(double(er),'fro') > realTol
                error('PhasorArray:sinCosForm:imagPartExceedsTol', 'Imaginary part norm (%e) exceeds tolerance (%e). Set isReal=false or increase tolerance.', norm(double(er),'fro'), realTol)
            end
            r=real(r);
        end
    end
    function r = AngleAmpForm(pA1)
        % ANGLEAMPFORM Convert PhasorArray to Angle-Amplitude form.
        %
        %   r = ANGLEAMPFORM(pA1) converts the PhasorArray `pA1` into
        %   the **Angle-Amplitude representation**, ensuring:
        %       - **Negative harmonics (`-k`) store phase angles (left side).**
        %       - **Positive harmonics (`+k`) store amplitudes (right side).**
        %       - The order `k` and `-k` remain **symmetric** around the DC component.
        %
        %   This format assumes `pA1` represents a **real-valued** periodic matrix.
        %
        %   Input:
        %       pA1 - (PhasorArray) The input PhasorArray object.
        %
        %   Output:
        %       r - (array) The Angle-Amplitude form representation:
        %           - Dimensions: `(m × n × (2h+1))`
        %           - Order `-k` maps to `(:,:,h-k+1)` (phase angle).
        %           - Order `+k` maps to `(:,:,h+k+1)` (amplitude).
        %
        %   Example:
        %       r = AngleAmpForm(pA1);
        %
        %   See also: SinCosForm, ImagRealForm
        h=pA1.h;
        r=cat(3,angle(flip(pA1.phas(1:h))),real(pA1.phas(0)),abs(pA1.phas(1:h)));
    end
    function r = CosSinForm(pA1)
        % COSSINFORM Convert PhasorArray to Cosine-Sine form.
        %
        %   r = COSSINFORM(pA1) converts the PhasorArray `pA1` into
        %   the **Cosine-Sine representation**, ensuring:
        %       - **Negative harmonics (`-k`) store cosine coefficients (left side).**
        %       - **Positive harmonics (`+k`) store sine coefficients (right side).**
        %
        %   Input:
        %       pA1 - (PhasorArray) The input PhasorArray object.
        %
        %   Output:
        %       r - (array) The Cosine-Sine form representation.
        %
        %   Example:
        %       r = CosSinForm(pA1);
        %
        %   See also: SinCosForm, ImagRealForm
        r = flip(SinCosForm(pA1),3);
    end

    function r = squeeval(pA1)% SQUEEVAL Evaluate and squeeze PhasorArray.
        %
        %   r = SQUEEVAL(pA1) evaluates the PhasorArray `pA1` and applies `squeeze`
        %   to remove singleton dimensions. This is useful for cases where `pA1`
        %   represents a scalar-valued PhasorArray.
        %
        %   Input:
        %       pA1 - (PhasorArray) The input PhasorArray object.
        %
        %   Output:
        %       r - (array) The squeezed numerical evaluation of `pA1.value`.
        %
        %   Example:
        %       r = squeeval(pA1);
        %
        %   See also: value, squeeze
        r=squeeze(value(pA1));
    end

    function r = expandBase(pA1,m)
        % EXPANDBASE Insert zeroed phasors to change the frequency base.
        %
        %   r = EXPANDBASE(pA1, m) modifies the PhasorArray `pA1` by inserting `m-1`
        %   zero-valued phasors between each existing phasor, effectively changing
        %   the frequency base from `w0` to `w1 = w0 / m`.
        %
        %   This transformation results in a PhasorArray with an increased harmonic
        %   resolution, suitable for frequency resampling or down-sampling applications.
        %
        %   Input:
        %       pA1 - (PhasorArray) The input PhasorArray.
        %       m  - (integer) The expansion factor (must be ≥ 1).
        %
        %   Output:
        %       r  - (PhasorArray) The modified PhasorArray with inserted zero phasors.
        %
        %   Behavior:
        %       - The harmonics are redistributed such that phasors appear only at
        %         indices `k*m`, while the in-between phasors are zeroed.
        %       - The resulting PhasorArray corresponds to the same periodic function,
        %         but expressed in a lower fundamental frequency `w1 = w0/m`.
        %
        %   Example:
        %       A_new = expandBase(A, 3); % Adjusts A for frequency w0/3
        %
        %   See also: squishBase, pad
        h = pA1.h;
        nx = size(pA1,1);
        ny = size(pA1,2);

        r = zeros(nx,ny,2*m*h+1);
        r(:,:,1:m:end) = pA1.value;
        r = PhasorArray(r);
    end

    function [r, normDeleted,normbase] = squishBase(pA1,m)
        % SQUISHBASE Remove phasors to change the frequency base.
        %
        %   r = SQUISHBASE(pA1, m) modifies the PhasorArray `pA1` by retaining only
        %   every `m`-th phasor, effectively changing the frequency base from `w0`
        %   to `w1 = w0 * m`.
        %
        %   This transformation results in a reduced harmonic resolution, increasing
        %   the fundamental frequency but potentially losing information in the process.
        %
        %   Input:
        %       pA1 - (PhasorArray) The input PhasorArray.
        %       m  - (integer) The reduction factor (must be ≥ 1).
        %
        %   Output:
        %       r  - (PhasorArray) The modified PhasorArray with retained harmonics.
        %
        %   Behavior:
        %       - Only harmonics at indices `k*m` are retained, while others are discarded.
        %       - If `pA1.h` is not a multiple of `m`, it is first padded with zeros.
        %       - A warning is issued if nonzero harmonics are removed, as this can
        %         lead to information loss and an inaccurate representation.
        %
        %   Example:
        %       A_new = squishBase(A, 3); % Adjusts A for frequency w0 * 3
        %
        %   See also: expandBase, pad
        h = pA1.h;
        nx = size(pA1,1);
        ny = size(pA1,2);
        %ensure pA1.h is a multiple of m, otherwise pad pA1 with zeros
        if mod(h,m)~=0
            pA1 = pad(pA1,m-mod(h,m));
            h = pA1.h;
        end
        normbase = norm(pA1(:),'fro');
        r = pA1(:,:,1:m:end);

        r = PhasorArray(r);

        %evaluate the norm of deleted phasors in original phasorArray
        deleted = pA1(:,:,setdiff(1:(2*h+1),1:m:(2*h+1)));
        normDeleted = norm(deleted(:),'fro');
        if normDeleted/normbase > 1e-4
            warning('PhasorArray:squishBase:truncationLoss', 'Deleted phasors have a Frobenius norm of %e (%.5f pct of total norm). The result may not accurately represent the original signal.', normDeleted, normDeleted/normbase*100)
        end
    end

end

methods (Access=protected)

    function header = getHeader(obj)
        nn = ndims(obj);

        %find class of obj.value to display in header
        classStr = sprintf('<a href="matlab:helpPopup %s">%s</a>',class(obj.value),class(obj.value));

        if nn<=3
            header = sprintf('%dx%dx%d <a href="matlab:helpPopup PhasorArray">PhasorArray</a> of %s',size(obj,1),size(obj,2),size(obj,3),classStr);
        else
            if nn<=5
                header = sprintf('%dx%d Array of %dx%dx%d <a href="matlab:helpPopup PhasorArray">PhasorArray</a> of %s',size(obj,4),size(obj,5),size(obj,1),size(obj,2),size(obj,3),classStr);
            else
                header = sprintf('%dD Array of %dx%dx%d <a href="matlab:helpPopup PhasorArray">PhasorArray</a> of %s',nn-3,size(obj,1),size(obj,2),size(obj,3),classStr);
            end
        end
        if nn>3
            endHeader = '';
        elseif isa(obj.value,'double') && iszero(obj)
            endHeader = sprintf(' representing a %dx%d <strong>zero-valued matrix</strong>',size(obj,1),size(obj,2));
        elseif isreal(obj)
            endHeader = sprintf(' representing a %dx%d <strong>real-valued</strong> periodic matrix with %d harmonics',size(obj,1),size(obj,2),obj.h);
        else
            endHeader = sprintf(' representing a %dx%d <strong>complex-valued</strong> periodic matrix with %d harmonics',size(obj,1),size(obj,2),obj.h);
        end

        header = sprintf('%s%s\n',header,endHeader);
    end

    function propgrp = getPropertyGroups(obj)
        propgrp = [];
    end
    function varargout = parenReference(obj,indexOp)
        obj.Phasor3D = obj.Phasor3D.(indexOp(1));
        if isscalar(indexOp)
            varargout{1} = pvalue(obj);
            return;
        end
        % Syntax for forwarding indexing operations
        [varargout{1:nargout}] = obj.(indexOp(2:end));
    end

    function obj = parenAssign(obj,indexOp,varargin)
        if isscalar(indexOp)
            assert(nargin==3);
            rhs = varargin{1};
            if ~isempty(obj)
                if isa(rhs,'PhasorArray')
                    obj.Phasor3D.(indexOp) = rhs.value;
                else
                    obj.Phasor3D.(indexOp) = rhs;
                end
            else
                %obj is uninitialized, assignement to undeclared variable eg A = B, A not initialized in workspace, B a phasorArray
                %so we want something like obj(indexOp) = rhs.value
                %and then obj = PhasorArray(obj)

                phs3D = rhs.value;

                %creacte a new array/whatever class nphs3D st nphs3D(indexOp) = phs3D
                %then obj = PhasorArray(nphs3D)
                nphs3D = zeros(size(phs3D));
                nphs3D(indexOp.Indices{:}) = phs3D;
                obj = PhasorArray(nphs3D);

            end
            return;
        end
        [obj.(indexOp(2:end))] = varargin{:};
    end

    function n = parenListLength(obj,indexOp,ctx)
        if numel(indexOp) <= 2
            n = 1;
            return;
        end
        containedObj = obj.(indexOp(1:2));
        n = listLength(containedObj,indexOp(3:end),ctx);
    end

    function obj = parenDelete(obj,indexOp)
        obj.Phasor3D.(indexOp) = [];
    end

    function  varargout=braceReference(obj,indexOp)
        if numel(indexOp(1).Indices)==3
            argou=obj.sub(indexOp(1).Indices{1:2});
            varargout{1}=argou.phas(indexOp(1).Indices{3});
        elseif numel(indexOp(1).Indices)<3
            [varargout{1:nargout}]=obj.sub(indexOp(1).Indices{:});
            if ~isscalar(indexOp)
                for ii=1:numel(varargout)
                    varargout{ii}=varargout{ii}.(indexOp(2:end));
                end
            end
        else
            argou = obj.sub(indexOp(1).Indices{1:2});
            argou = value(argou.phas(indexOp(1).Indices{3}));
            ind = repmat({':'},1,ndims(argou));
            [ind{4:end}] = indexOp(1).Indices{4:end};
            varargout{1} = argou(ind{:});
            if strcmp(indexOp(1).Indices{3},':')
                varargout{1} = PhasorArray(varargout{1});
            end
        end
    end

    function n = braceListLength(obj,indexOp,indexContext)
        if numel(indexOp) <= 2
            n = 1;
            return;
        end
        n=1;
    end

    function obj = braceAssign(obj,indexOp,varargin)

        h1=size(obj,3);
        h2=size(varargin{:},3);
        m1=size(obj,1);
        m2=size(obj,2);

        if numel(indexOp.Indices)==3
            n1=indexOp.Indices{1};
            n2=indexOp.Indices{2};

            n3=indexOp.Indices{3};
            if ischar(n3) % %we get ":" , so necessarly a phasorArray
                obj(n1,n2,n3)=varargin{:};
                return
            end
            if max(abs(n3))>obj.h
                obj=obj.pad(max(abs(n3))-obj.h);
            end
            n3=n3+1+obj.h;
            % Writing k alone leaves -k untouched, which silently costs the
            % realness of A(t). Checked on the written pair only: global
            % isreal is O(n*m*h) and would double the cost of an assignment.
            mir  = 2*obj.h + 2 - n3;
            wasC = pairIsConjugate(obj, n3, mir);
            obj(n1,n2,n3)=varargin{:};
            if wasC && ~pairIsConjugate(obj, n3, mir)
                warning('PhasorArray:braceAssign:brokenConjugatePair', ...
                    ['Harmonics %s and their opposites are no longer conjugate, so A(t) ' ...
                     'is no longer real. Assign the mirror too, or use setHarmonic.'], ...
                    mat2str(n3 - 1 - obj.h));
            end
            return
        end

        if h1<h2
            obj=obj.pad((h2-h1)/2);
        else
            varargin{:}=PhasorArray(PhasorArrayPad(varargin{:},(h1-h2)/2));
        end

        if isscalar(indexOp)
            switch numel(indexOp.Indices)
                case 1
                    m=indexOp.Indices{1};
                    if  numel(m)==obj.numelt %numelt out numel (A(t)) is size(A,1) times size(A,2) returns the dimension of the temporal matrix i.e. n1 \times n2, so it would be a logical array as input
                        if any((m ~= 0) & (m ~= 1)) % check m is logical
                            error('PhasorArray:subsref:invalidIndex', 'Array indices must be positive integers or logical values.')
                        end

                        n1=find(indexOp.Indices{1});
                    else
                        n1=indexOp.Indices{1};
                    end
                    n2=1;
                    % obj
                    obj=obj{:};
                case 2
                    n1=indexOp.Indices{1};
                    n2=indexOp.Indices{2};
            end
            if isscalar(varargin{:})
                varargin{:}=varargin{:}*ones(numel(n1),numel(n2));
            end

            %                 obj(n1,n2,:) = vect(varargin{:});
            obj(n1,n2,:) = (varargin{:});
            obj=reshape(obj,m1,m2,[]);
            return;
        end
        %             [braceReference(obj,indexOp)] = varargin{:};
    end

end

methods (Access=public)

    function ind = end(obj,k,n)
        s1 = size(obj,1);
        s2 = size(obj,2);
        s3 = size(obj,3);
        sz = [s1 s2 s3];
        if n == 1
            ind = s1*s2;
            return
        end
        ind = sz(k);
    end

    function out = value(obj)
        %VALUE Raw coefficient array, [n x m x (2h+1)], harmonics ascending.
        %
        %   Slice k holds harmonic k-h-1, so the centre slice is the DC component.
        %
        %   See also: pvalue, phas.
        out = obj.Phasor3D;
    end

    function out = sum(obj,dim)
        %SUM Compute the sum along a given dimension.
        %
        %   out = SUM(obj, dim) returns the sum along the specified dimension.
        %   If dim is 1 or 2, the output remains a PhasorArray.
        arguments
            obj
            dim=1
        end
        out = sum(obj.Phasor3D,dim);
        if dim==1 || dim==2
            out=PhasorArray(out);
        end
    end

    function out = cat(dim,varargin)
        %CAT Concatenate multiple PhasorArrays along a given dimension.
        %
        %   out = CAT(dim, A1, A2, ...) concatenates multiple PhasorArrays or
        %   numeric arrays along the specified dimension.
        numCatArrays = nargin-1;
        if dim==1 || dim==2
            varargin2=PhasorUnif(varargin{:});
        else
            varargin2=varargin;
        end

        newArgs = cell(numCatArrays,1);
        for ix = 1:numCatArrays
            if isa(varargin2{ix},'PhasorArray')
                newArgs{ix} = varargin2{ix}.Phasor3D;
            else
                newArgs{ix} = varargin2{ix};
            end
        end
        out = PhasorArray(cat(dim,newArgs{:}));
    end

    function varargout = size(obj,varargin)
        [varargout{1:nargout}] = size(obj.Phasor3D,varargin{:});
    end

    function r = sdpval(pA1)
        %SDPVAL Extract the numerical value of a PhasorArray whose value is an SDPVAR.
        %
        %    r = sdpval(pA1) extract the numerical value of pA1 a phasorArray with NDSDPVAR 3DArray.
        %    Basically perform PhasorArray(value(value(pA1)))
        if  isspecial(pA1)
            pA1=pA1.value;
        end
        r=PhasorArray(value(pA1));
    end

    [K,P,res,info] = place(A,B,poles,nvp) ;


end

methods
    function r = Value(obj)
        %VALUE Alias of value, kept for the call sites that use the capitalised name.
        %
        %   Prefer value. MATLAB is case-sensitive, so the two names are genuinely
        %   distinct methods and appear separately in help.
        %
        %   See also: value.
        r=obj.value;
    end
end

methods (Static, Access=public)
    % Implementation in @PhasorArray/randomWithNPole.m
    [Aper, info] = randomWithNPole(J_or_V, h, varargin)
    % Implementation in @PhasorArray/randomSPD.m
    P = randomSPD(nx, h, varargin)
    % Implementation in @PhasorArray/randomPhasorArrayWithPole.m
    obj = randomPhasorArrayWithPole(nx, poles, T, varargin)
    % Implementation in @PhasorArray/sym.m
    A = sym(nx, ny, h, name, varargin)
    % Implementation in @PhasorArray/ndsdpvar.m
    [P, held] = ndsdpvar(n1, n2, h, varargin)

    function out = builtin(funcName, varargin)
        % Intercepts function-style calls: methodName(obj)
        out = builtin(funcName, varargin{:});
    end

    function obj = scalar(varargin)
        % SCALAR Create a scalar PhasorArray.
        %
        %   obj = SCALAR(varargin) creates a scalar PhasorArray from the input arguments.
        %   The input can be:
        %       - A scalar value.
        %       - A vector of phasor coefficients.
        %       - A 3D array representing the harmonics.
        %
        %   Inputs:
        %       varargin:
        %           - 1 argument: A vector, which will be converted into 3D array.
        %           - 2 arguments: The first input is a scalar (DC component of phasor array),
        %             and the second is a vector, the positive phasors of the signal.
        %               The optional argument `isreal` is set to `true` by default if not provided.
        %       nvp:
        %           - reduce: Boolean flag indicating whether to reduce the array (default: false).
        %           - isreal: Boolean flag to enforce positive parts of the phasor (default: false).
        %
        %   Outputs:
        %       obj: A PhasorArray object containing the scalar phasor data.
        %
        %   Example:
        %       obj = PhasorArray.scalar(5); % Creates a scalar PhasorArray with a DC component of 5.
        %       obj = PhasorArray.scalar([1, 2, 3]); % Creates a scalar PhasorArray with specified phasors.
        %       obj = PhasorArray.scalar(rand(1, 1, 5)); % Creates a scalar PhasorArray from a 3D array.
        %
        %   See also: PhasorArray, ScalarPhasorArray

        obj = ScalarPhasorArray(varargin{:});
    end

    function obj = fromTBMatrix(TBMatrix, nx, size_type)
        % FROMTBMATRIX Create a PhasorArray from a Toeplitz Block matrix.
        %   obj = FROMTBMATRIX(TBMatrix, nx, size_type) converts a Toeplitz Block matrix
        %   into a PhasorArray by extracting the anti-diagonal elements.
        %
        %   Inputs:
        %     TBMatrix  - (matrix) The Toeplitz Block matrix to be converted.
        %     nx        - (scalar) The first dimension of the time-dependent matrix A(t).
        %     size_type - (string) Specifies which dimension of A(t) is provided.
        %              - 'n1': The first dimension is nx.
        %              - 'n2': The second dimension is nx.
        %               Default: 'n1'.
        %
        %   Outputs:
        %     obj - (PhasorArray) The converted PhasorArray.
        %
        %   Example:
        %     TBMatrix = rand(15, 15); % Example Toeplitz Block matrix
        %     phasorArray = PhasorArray.fromTBMatrix(TBMatrix, 5, 'n1');
        %       since 5 is provided, we deduce here that each block is of size 3x3, and second dim of A(t) is 3.
        %       so truncature is done to (3-1)/2=1;
        %
        %   See also: TB2array, TB
        pA = TB2array(TBMatrix, nx, size_type);
        obj = PhasorArray(pA);
    end

    function obj = fromF_tb(F_tb, n1, n2,nvp)
        % FROMF_tb Create a PhasorArray from a Fourier Transform (TB format).
        %   obj = FROMF_tb(F_tb, n1, n2) converts a vector representing the Fourier
        %   decomposition in a Toeplitz block (TB) format into a PhasorArray.
        %
        %   Inputs:
        %     F_tb - (vector) A vector containing the Fourier coefficients in TB format.
        %     n1    - (scalar) The first dimension of the matrix.
        %     n2    - (scalar,option) The second dimension of the
        %     matrix, argument n2 is to be used to convert
        %     F_tb(col(M)) back to phasor array M, otherwise
        %
        %   Outputs:
        %     obj - (PhasorArray) The converted PhasorArray.
        %
        %   Example:
        %     M = PhasorArray.random(3,4, 8);
        %    TF_M = F_tb(M,5)% is a 4*(2*5+1) times 4 matrix
        %      M_r = PhasorArray.fromF_tb(TF_M,3) or PhasorArray.fromF_tb(TF_M,"h",5)
        %
        %   TF_colM = F_tb(M{:},12) is a (4*3*(2*12+1)) times 1 col matrix
        %      M_r_fromCol = PhasorArray.fromF_tb(TF_colM,3,4)
        %
        %   See also: F_tb_2_PhasorArray,F_tb
        arguments
            F_tb
            n1 =[]
            n2 =[]
            nvp.h =[]
        end
        if isempty(n2)
            n2 = size(F_tb,2);
        end

        if ~isempty(nvp.h)
            n1h = size(F_tb,1)/(2*nvp.h+1);
            if ~isempty(n1) && (n1~=n1h || (n1*n2)~=n1h)
                error('PhasorArray:fromF_tb:incompatibleArguments', 'n1=%d, n2=%d, h=%d are not compatible with the input dimensions.', n1, n2, nvp.h)
            end
        end
        PA = F_tb_2_PhasorArray(F_tb, n1, n2);
        obj = PhasorArray(PA);
    end

    function out = time2Phasor(At,nT,t,nvp)
        %TIME2PHASOR Convert a time-dependent matrix into a PhasorArray representation.
        %
        %   out = TIME2PHASOR(At, nT, t, nvp) converts a 3D array representing a
        %   real, time-dependent matrix (with time stored along the 3rd dimension)
        %   into a PhasorArray by computing its Fourier coefficients.
        %
        %   Assumptions:
        %     - The sampling is uniform: Ts = t(2) - t(1).
        %     - The final time matches a full number of periods: t(end) = nT*T - dt.
        %
        %   Inputs:
        %     At       - (3D array) Time-dependent matrix, stored as MxNxtime.
        %     nT       - (integer, default: 1) Number of periods captured in `t`.
        %     t        - (vector, optional) Time vector. If empty, it is inferred.
        %     nvp     - (Optional) Name-value pair arguments:
        %         'truncIndex' (integer, default: Inf)  - Maximum harmonic order for truncation.
        %         'real'       (logical, default: true) - Enforce conjugate symmetry (A_k = conj(A_-k)).
        %         'timeDim'    (integer, default: 3)    - Dimension along which time varies.
        %
        %   Outputs:
        %     out - (PhasorArray) The computed PhasorArray representation of `At`.
        %
        %   Example:
        %     t = linspace(0, 2*pi, 100);
        %     At = cat(3, sin(t), cos(t));
        %     Ph = time2Phasor(At, 1, t, 'truncIndex', 5);
        %
        %   See also: PhasorArray, TimeArray2Phasors
        %
        arguments
            At
            nT=1
            t=[]
            nvp.truncIndex=Inf
            nvp.real=true
            nvp.timeDim=3
        end
        out = PhasorArray(TimeArray2Phasors(At,nT,t,'truncIndex',nvp.truncIndex,'isReal',nvp.real,'timeDim',nvp.timeDim));
    end

    function [phasorArr,feval] = funcToPhasorArray(func, T, n,nvp)
        %FUNCTOPHASORARRAY Convert a time function into a PhasorArray.
        %
        %   [phasorArr, feval] = FUNCTOPHASORARRAY(func, T, n, nvp) evaluates a
        %   time-dependent function over a uniform time grid and converts it into
        %   a PhasorArray using a Fourier transform.
        %
        %   Inputs:
        %     func - (function_handle) Time-dependent function, f(t), returning an MxN matrix.
        %     T    - (scalar, default: 1) Period of the function.
        %     n    - (integer, default: 4) Number of frequency bins used (2^n time steps).
        %     nvp - (Optional) Name-value pair arguments:
        %         'reduce'     (logical, default: true)  - Reduce the output PhasorArray.
        %         'reduceTol'  (double, default: 1e-15) - Threshold for reducing small harmonics.
        %         'isReal'     (logical, default: isreal(func(t))) - Enforce the
        %                      conjugate symmetry Ak = conj(A-k). Set it false for
        %                      a complex-valued func, true to clean the round-off
        %                      asymmetry of a real one.
        %
        %   Outputs:
        %     phasorArr - (PhasorArray) The computed PhasorArray representation of `func(t)`.
        %     feval     - (struct) Structure containing function evaluations for debugging:
        %         - feval.T  : Period used.
        %         - feval.n  : Number of frequency bins.
        %         - feval.dt : Time step.
        %         - feval.func : Function handle used.
        %         - feval.At  : Evaluated matrix over time.
        %         - feval.t   : Time vector used.
        %
        %   Example:
        %     func = @(t) [sin(t); cos(t)];
        %     T = 2*pi;
        %     n = 5;
        %     [Ph, fData] = funcToPhasorArray(func, T, n);
        %
        %   See also: time2Phasor, PhasorArray
        %
        arguments
            func (1,1) {mustBeA(func, 'function_handle')} = @(t) ones(3,3)
            T (1,1) {mustBeNumeric, mustBePositive} = 1
            n (1,1) {mustBeNumeric, mustBePositive, mustBeInteger} = 4
            nvp.reduce = true;
            nvp.reduceTol = 1e-15;
            nvp.isReal = [];
        end

        % Calculate dt
        dt = T / (2^n);

        % Create a time vector
        t = 0:dt:T-dt;

        % Evaluate the function over the time vector
        At = arrayfun(func, t, 'UniformOutput', false);
        At = cat(3, At{:});

        AtT = func( T);
        
        %Redacted : function can be periodic and discontinuous, so first
        %test is meaningleass

        %evaluate the jump between the last and first value
        % if norm(At(:,:,end)-At(:,:,1),'fro')>1e-10
        %     warning('The function has discontinuities, or a steep derivative at the end, or has a jump between the last and first value, the result may be incorrect. Increase sampling (m value). Jump value is %d (froebenius norm of f(T-dt)-f(0))',norm(At(:,:,end)-At(:,:,1),'fro'));
        % end
        % 
        if norm(At(:,:,1)-AtT,'fro')>1e-10
            warning('PhasorArray:fromFun:notPeriodic', 'A(T) ~= A(0): jump Frobenius norm = %e. The function may not be T-periodic; results may be inaccurate.', norm(At(:,:,1)-AtT,'fro'));
        end

        % time2Phasor enforces Ak = conj(A-k) by default, which silently
        % destroys a complex-valued func, so the flag follows the samples.
        if isempty(nvp.isReal)
            nvp.isReal = isreal(At);
        end
        phasorArr = PhasorArray.time2Phasor(At, 'real', nvp.isReal);
        if nvp.reduce
            phasorArr = reduce(phasorArr, "reduceThreshold",nvp.reduceTol,"reduceMethod","relative","exclude0Phasor",false,"hardThresholdPhasors",true);

        end
        if nargout>1
            feval=struct;
            feval.T=T;
            feval.n=n;
            feval.dt=dt;
            feval.func=func;
            feval.At=At;
            feval.t=t;
            feval.plot = @(n,alone) At2plot(n,alone);
        end



        function At2plot(n,alone)
            %PLOTFUNC Plot time-domain representation of function evaluations.
            %
            %   plotFunc(numPeriods, mode) plots the function evaluations over `numPeriods`
            %   periods. The `mode` argument can be 'alone' (plot time-domain only) or 'both'
            %   (overlay with the PhasorArray reconstruction).
            %
            %   Example:
            %     feval.plot(4, 'both');
            %
            if nargin==0
                n=4;
                alone="alone";
            end

            if nargin==1
                if (ischar(n) || isstring(n))
                    alone=n;
                    n=4;


                elseif nargin==1 && isnumeric(n)
                    alone="alone";
                end
            end

            %assert if alone is "alone" or "both"
            assert(strcmp(alone,"alone")||strcmp(alone,"both"),"alone must be 'alone' or 'both'");
            dt=T/(100);
            t = 0:dt:n*T;
            At = arrayfun(func, t, 'UniformOutput', false);
            At = cat(3, At{:});

            %plot each A(i,j,t) as a function of time
            for i=1:size(At,1)
                for j=1:size(At,2)
                    subplot(size(At,1),size(At,2),(i-1)*size(At,2)+j);
                    plot(t,squeeze(At(i,j,:)));
                    title(sprintf('A_{%d,%d}(t)',i,j));
                end
            end
            if strcmp(alone,"both")
                hold on;
                phasorArr.plot(T,[0 n*T]);
                hold off;
            end
        end
    end


    function out = cqt2ScalarPhasor(cqtobj,nvp)
        %CQT2SCALARPHASOR Convert a CQT object to a ScalarPhasorArray representation.
        %
        %   out = CQT2SCALARPHASOR(cqtobj, nvp) extracts the symbolic Toeplitz
        %   representation of a compact quasi-Toeplitz (CQT) object and maps it
        %   to a ScalarPhasorArray.
        %
        %   Inputs:
        %     cqtobj - (CQT object) The input compact quasi-Toeplitz object.
        %     nvp   - (Optional) Name-value pair arguments:
        %         'isReal' (logical, default: false) - If true, enforces real-valued output
        %                                              by ensuring conjugate symmetry in
        %                                              the phasor representation.
        %
        %   Outputs:
        %     out - (ScalarPhasorArray) The extracted and converted phasor representation.
        %
        %   Example:
        %     cqtA = cqt(rand(10), rand(10)); % Example CQT matrix
        %     phasorA = cqt2ScalarPhasor(cqtA, 'isReal', true);
        %
        %   See also: ScalarPhasorArray, cqt
        %
        arguments
            cqtobj %object of class cqt
            nvp.isReal = false %enforce real output (ie conjugate phasore pos/neg)
        end
        [nsTA,psTA]=symbol(cqtobj);
        n=min(numel(nsTA),numel(psTA));
        nsTA=nsTA(1:n);
        psTA=psTA(1:n);
        out=ScalarPhasorArray([flip(nsTA(2:end)) psTA]);
        if nvp.isReal
            out = mreal(out);
        end
    end


    function obj = empty(varargin)
        %EMPTY Create an empty PhasorArray of specified size.
        %   See also: zeros, ones, eye
        obj = PhasorArray(double.empty(varargin{:}));
    end
    function obj = zeros(varargin)
        %ZEROS Create a zero-filled PhasorArray of specified size.
        %   See also: ones, eye, empty
        obj = PhasorArray(zeros(varargin{:}));
    end
    function obj = ones(varargin)
        %ONES Create a PhasorArray filled with ones.
        %   See also: zeros, eye, empty
        obj = PhasorArray(ones(varargin{:}));
    end
    function obj = eye(varargin)
        %EYE Create a PhasorArray with identity matrices at all harmonics.
        %   EYE(n) returns an `n × n × 1` PhasorArray filled with identity matrices.
        %   EYE(n, m) returns an `n × m × 1` PhasorArray filled with identity matrices.
        %   EYE(n, m, h) returns an `n × m × (2h+1)` PhasorArray where each slice along the 3rd dimension is an identity matrix.
        %   Unlike a standard identity matrix, this function replicates the identity structure across all harmonic indices.
        %   See also: zeros, ones, empty
        if nargin==1
            varargin{2}=varargin{1};
            varargin{3}=0;
        elseif nargin==2
            varargin{3}=0;
        elseif nargin==3

        else
        end
        u=repmat(eye(varargin{1:2}),[1 1 2*varargin{3}+1]);
        obj = PhasorArray(u);
    end
    function obj = random(nx,ny,h,nvp)
        % RANDOM Generate a random PhasorArray with optional structure constraints.
        %
        %   OBJ = RANDOM(NX, NY, H, ARG) generates a random 3D PhasorArray representing
        %   a time-periodic matrix with harmonics and structural constraints.
        %
        %   Inputs:
        %     NX  - (integer) Number of rows.
        %     NY  - (integer) Number of columns.
        %     H   - (integer) Number of harmonics.
        %     ARG - (optional) Name-value pair arguments:
        %             - 'symmetry' (string array) Symmetry class, any combination of
        %                   the 14 names listed in PHASORSYMMETRY (default: "real").
        %                   Pass [] for a full complex array.
        %             - 'T'           (scalar) Time period (default: `2*pi`).
        %             - 'average_power_decay' (scalar) Power decay rate of harmonics (default: `2`).
        %
        %   Outputs:
        %     OBJ - (PhasorArray) Randomly generated PhasorArray.
        %
        %   'symmetry' constrains structure, never the spectrum. For prescribed Floquet
        %   exponents use randomWithNPole; for A(t) positive definite, randomSPD.
        %
        %   Example
        %       A = PhasorArray.random(3, 3, 5, "symmetry", ["real" "symmetric"]);
        %       B = PhasorArray.random(3, 3, 5, "symmetry", []);      % full complex
        %
        %   See also: rand_phasor, phasorSymmetry, PhasorArray.randomSPD,
        %             randomPhasorArrayWithPole
        arguments
            nx
            ny
            h
            nvp.symmetry string = "real"
            nvp.T=2*pi
            nvp.average_power_decay=2;
        end
        if nargin ==1
            switch numel(nx)
                case 1
                case 2
                    ny = nx(2);
                    nx = nx(1);
                    h = 0;
                case 3
                    ny = nx(2);
                    h  = nx(3);
                    nx = nx(1);
            end
        end
        nvp.output ='PhasorArray'; %or PhasorArray

        C=namedargs2cell(nvp);
        obj = rand_phasor(nx,ny,h,C{:});
    end
    function out = cos(phas,order)
        % COS Produce PhasorArray representation of cosine function.
        %
        %   out = COS(phas, order) returns a PhasorArray representation of the
        %   cosine function dephased by `phas` and with a specified harmonic `order`.
        %
        %   Inputs:
        %     phas  - (scalar, optional) Phase shift of the cosine function. Default is 0.
        %     order - (integer, optional) Harmonic order of the cosine function. Default is 1.
        %
        %   Outputs:
        %     out - (PhasorArray) PhasorArray representation of the cosine function.
        %
        %   Example:
        %     out = PhasorArray.cos(); % Default cosine function with no phase shift and order 1, ie cos(ω t).
        %     out = PhasorArray.cos(pi/4, 2); % Cosine function with phase shift of pi/4 and harmonic order 2, ie cos(2ω t + pi/4).
        %     out = PhasorArray.cos([0 -2π/3 2π/3]) % Cosine function with phase shift of 0, -2π/3, and 2π/3, ie [cos(ω t), cos(ω t - 2π/3), cos(ω t + 2π/3)].
        %
        %   See also: sin
        arguments
            phas = 0
            order {mustBePositive, mustBeInteger} = 1
        end
        assert(order>0,'order of the cosine must be positive')
        n0 = zeros(1,2*order-1);
        out = PhasorArray.scalar([1/2 n0 1/2]);
        out = out.PhaseShift(phas);
    end

    function out = sin(phas,order)
        % SIN Produce PhasorArray representation of sine function.
        %
        %   out = SIN(phas, order) returns a PhasorArray representation of the
        %   sine function dephased by `phas` and with a specified harmonic `order`.
        %
        %   Inputs:
        %     phas  - (scalar, optional) Phase shift of the sine function. Default is 0.
        %     order - (integer, optional) Harmonic order of the sine function. Default is 1.
        %
        %   Outputs:
        %     out - (PhasorArray) PhasorArray representation of the sine function.
        %
        %   Example:
        %     out = PhasorArray.sin(); % Default sine function with no phase shift and order 1, ie sin(ω t).
        %     out = PhasorArray.sin(pi/4, 2); % Sine function with phase shift of pi/4 and harmonic order 2, ie sin(2ω t + pi/4).
        %     out = PhasorArray.sin([0 -2π/3 2π/3]) % Sine function with phase shift of 0, -2π/3, and 2π/3, ie [sin(ω t), sin(ω t - 2π/3), sin(ω t + 2π/3)].
        %
        %   See also: cos
        arguments
            phas = 0
            order {mustBePositive, mustBeInteger} = 1
        end
        assert(order>0,'order of the sine must be positive')
        n0 = zeros(1,2*order-1);
        out = PhasorArray.scalar([-1/(2*1i) n0 1/(2*1i)]);
        out = out.PhaseShift(phas);
    end

    % Coordinate transforms — implementations in @PhasorArray/*.m
    out = Rotdq0(dephase, order, include0)
    out = Concordia(include0)
    out = Clark(include0)
    out = Park(dephase, order, include0)
    out = negativeRotdq0(dephase, order, include0)
    out = negativePark(dephase, order, include0)
    out = zeroPosNegSequenceDQ(dephase, order, include0)
    out = dq0(dephase, order, include0)
    out = negativeDQ0(dephase, order, include0)
    out = ZPNSequence(dephase, order, include0)

end
end

%% =========================================================================
function [Eew, E] = energyOf(X, elementwise, nOut)
%ENERGYOF  Squared l2 norm of a PhasorArray over its harmonic dimension.
%
%   SYNTAX
%     [Eew, E] = energyOf(X, elementwise, nOut)
%
%   DESCRIPTION
%   Shared body of realEnergy / imagEnergy / energy / ACenergy / DCenergy,
%   which differ only in which projection of the array is passed in.
%
%     E(i,j) = sum_k |X(i,j,k)|^2         (element-wise)
%     E      = sum_ij E(i,j)              (total)
%
%   By Parseval this equals the time-domain energy int |X(t)|^2 dt over one
%   period, so the same kernel serves every member of the family.
%
%   INPUTS
%     X           - PhasorArray (already projected by the calling method).
%     elementwise - Logical. With a single output, return the element-wise
%                   matrix instead of the total.
%     nOut        - The caller's nargout. With two outputs requested, Eew is
%                   always element-wise and E always the total.
%
%   OUTPUTS
%     Eew - Element-wise energy, or the total when nOut < 2 and ~elementwise.
%     E   - Total energy.
Eew = sum(X.value .* conj(X.value), 3);
E   = sum(Eew, 'all');
if nOut < 2 && ~elementwise
    Eew = E;
end
end

%% =========================================================================
function [r, R] = projectionMatch(pA1, projFcn, symSign, tol, label)
%PROJECTIONMATCH  Shared body of isreal (symSign=+1) and isimag (symSign=-1).
%
%   FLOAT FAST PATH. The question is a symmetry of the harmonic coefficients:
%       A(t) real           <=>  A_{-k} =  conj(A_k)
%       A(t) purely imag.   <=>  A_{-k} = -conj(A_k)
%   Fast path avoids building mreal/mimag. The tolerance is RELATIVE.
%   sym/sdpvar keep the projection route: they cannot go through max/abs.

if ~isspecial(pA1)
    v     = pA1.value;
    scale = max(abs(v), [], 'all');
    if isempty(scale) || scale == 0
        R = true(size(v));  r = true;  return
    end
    if isempty(tol), tol = 1e-12; end          % ismembertol's double default
    R = abs(v - symSign * conj(flip(v, 3))) <= tol * scale;
    r = all(R, 'all');
    return
end

undetermined = false;
try
    o1_p = projFcn(pA1);
catch
    undetermined = true;   % projection itself failed (was an unguarded throw)
end

% Only sym/sdpvar reach here: the float case returned from the fast path above.
if ~undetermined
    if isempty(tol), tol = eps; end
    try
        r1 = abs(real(pA1.value) - real(o1_p.value)) < tol;
        r2 = abs(imag(pA1.value) - imag(o1_p.value)) < tol;
    catch
        undetermined = true;
    end

    if ~undetermined && isa(r1, 'sym')
        try
            r1 = logical(r1);
            r2 = logical(r2);
        catch
            % Truth value depends on assumptions: ask isAlways instead.
            ws = warning('off', 'symbolic:sym:isAlways:TruthUnknown');
            try
                r1 = isAlways(r1);
                r2 = isAlways(r2);
            catch
                undetermined = true;
            end
            warning(ws);
        end
    end
end

if undetermined
    warning(['PhasorArray:' label ':undetermined'], ...
        ['%s could not be evaluated for a %s-backed PhasorArray; returning false. ', ...
         'This means "cannot determine", not "no".'], label, class(pA1.value));
    R = false(size(pA1.value));
    r = 0;
    return
end

R = zeros(size(pA1.value), 'logical');
R((r1 + r2) == 2) = true;
r = all(R, 'all');
end

%% =========================================================================
function tf = pairIsConjugate(obj, kIdx, mIdx)
%PAIRISCONJUGATE  Do slices kIdx and mIdx hold conjugate coefficients?
%   Local test, O(n*m), used to spot a brace assignment that breaks the
%   realness of A(t). Out-of-range mirrors read as broken.
if any(mIdx < 1) || any(mIdx > size(obj, 3))
    tf = false; return
end
v  = pvalue(obj);
tf = ~any(abs(v(:,:,kIdx) - conj(v(:,:,mIdx))) > 1e-12, 'all');
end

%% =========================================================================
function m = mollifierFT(xi)
%MOLLIFIERFT  Fourier transform of the unit bump, at the given frequencies.
%   phi(t) = exp(-1/(1-t^2)) on |t|<1, normalised. Real and even, so the
%   transform reduces to a cosine integral. All derivatives of phi vanish at
%   the endpoints, so the midpoint rule converges spectrally: 4096 points
%   already agree with 2^20 to nine digits.
persistent t phi Z
if isempty(t)
    N   = 4096;
    t   = -1 + (0.5:N-0.5)*(2/N);
    phi = exp(-1./(1-t.^2));
    Z   = sum(phi);
end
m = reshape(cos(xi(:) * t) * phi(:) / Z, size(xi));
end

%% =========================================================================
function [r, frac, crenel, Cph] = compareInTime(pA1, pA2, op)
%COMPAREINTIME  Order two periodic matrices pointwise in time.
%   r(i,j) holds at every sampled t; frac(i,j) is the fraction of the period
%   over which it holds, so r == (frac == 1). crenel is the indicator,
%   [n x m x Nt] logical on theta = (0:Nt-1)/Nt*2*pi; Cph is the same as a
%   PhasorArray of order Nt/2-1. MINUS supplies padding and broadcasting.
%
%   The difference must be real: MATLAB's ordering compares real parts on
%   complex input and drops the rest silently, so a complex one is refused.
%
%   Cph is a square wave. Its coefficients decay as 1/k and 46% of the energy
%   lies outside DC at any grid size, so only its DC term -- which is frac --
%   survives truncation.
%
%   The difference must be real in time: MATLAB's ordering operators compare
%   real parts on complex input and drop the imaginary half without warning,
%   which is the defect this function exists to remove. A complex difference
%   has no order, so it is refused rather than answered.
%
%   The verdict is sampled. D(t) is band-limited at h, so 2h+1 points already
%   determine it, but its extrema over t are not determined by samples alone;
%   the grid is oversampled well past that to make a near-boundary case
%   unlikely rather than impossible. Use a strict margin if the distinction
%   matters, or evaluate with evalp on your own grid.
D = pA1 - pA2;
% Ordering samples A(t)-B(t) on a grid, so it needs numbers. A decision
% variable has no order until it is solved, and a symbolic one has none at all.
if ~isnumeric(pvalue(D))
    error('PhasorArray:compare:symbolicPayload', ...
        ['Ordering compares sampled values and needs a numeric array; this ' ...
         'difference holds %s. Solve first and compare sdpval of the result, ' ...
         'or state the relation as a constraint rather than a comparison.'], ...
        class(pvalue(D)));
end
if ~isreal(D)
    error('PhasorArray:compare:complexOrder', ...
        ['Ordering requires A(t)-B(t) to be real; this difference is complex. ' ...
         'Compare mreal(A) with mreal(B), or use == / iszero for equality.']);
end
Nt = 2^nextpow2(max(256, 16*(D.h+1)));
th = (0:Nt-1)/Nt * 2*pi;
C  = op(real(evalp(D, th)), 0);
r  = all(C, 3);
if nargout > 1, frac   = mean(double(C), 3);                        end
if nargout > 2, crenel = C;                                         end
if nargout > 3, Cph    = PhasorArray(TimeArray2Phasors(double(C), 1, th)); end
end

%% =========================================================================
function r = pageInfNorm(X)
%PAGEINFNORM  norm(.,inf) page by page, without pagenorm.
%   pagenorm is R2022b while the class promises R2021b, so iszero used to fail
%   on the very releases the header advertises. norm(.,inf) switches semantics:
%   max(abs) on a vector page, largest absolute row sum on a matrix page.
%   Verified against pagenorm to 4.4e-16 over 500 random shapes.
if size(X,1) == 1 || size(X,2) == 1
    r = max(abs(X), [], [1 2]);
else
    r = max(sum(abs(X), 2), [], 1);
end
end

%% =========================================================================
function tf = symbolicSliceIsZero(S, tol)
%SYMBOLICSLICEISZERO  True if every entry of a symbolic slice is within tol of 0.
%   Mirrors the max-abs semantics of pagenorm(...,"inf") <= tol without
%   requiring a floating-point array. Returns false when the comparison cannot
%   be reduced to a truth value.
try
    c = abs(S) <= tol;
    if isa(c, 'sym')
        try
            c = logical(c);
        catch
            ws = warning('off', 'symbolic:sym:isAlways:TruthUnknown');
            c  = isAlways(c);
            warning(ws);
        end
    end
    tf = all(c(:));
catch
    tf = false;
end
end
