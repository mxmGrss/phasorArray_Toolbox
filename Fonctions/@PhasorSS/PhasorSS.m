% Class definition for PhasorSS, a class representing a state-space model with phasor arrays.
%
% Properties:
%   A - State matrix (PhasorArray)
%   B - Input matrix (PhasorArray)
%   C - Output matrix (PhasorArray)
%   D - Feedthrough matrix (PhasorArray)
%   p - LPV parameter for phase used for periodicity
%   isLPV - Logical flag indicating if the system is LPV or LTV
%   T - Period of the system if set to LTV
%   isReal - Logical flag indicating if the system is real
%   StateName - Cell array of state names
%   StateUnit - Cell array of state units
%   InputName - Cell array of input names
%   InputUnit - Cell array of input units
%   InputGroup - Structure of input groups
%   OutputName - Cell array of output names
%   OutputUnit - Cell array of output units
%   OutputGroup - Structure of output groups
%   Name - Name of the system
%   Notes - Cell array of notes
%   UserData - User-defined data
%   providedD - Hidden property to store the providedD and manage the empty D matrix case
%   providedC - Hidden property to store the providedC and manage the empty C matrix case (default is identity matrix)
%
% Methods:
%   % Info methods
%   size - Size of the PhasorSS object
%   isStatic - Check if the PhasorSS object is static (single D matrix, periodic)
%   nx - Number of states of the PhasorSS object
%   nu - Number of inputs of the PhasorSS object
%   ny - Number of outputs of the PhasorSS object
%
%   % Simulation methods
%   initial - Perform initial condition simulation of the PhasorSS object
%   step - Perform step input simulation of the PhasorSS object
%   impulse - Perform impulse input simulation of the PhasorSS object
%   lsim - Perform simulation of the PhasorSS object, with input history (u,t) and initial condition x0
%   lsimplot - Perform simulation of the PhasorSS object, with input history (u,t) and initial condition x0
%   stepu - Perform step input simulation of the PhasorSS object with a periodic input u
%
%   % Edition methods
%   setA - Set the A matrix of the PhasorSS object and update the isReal property if needed
%   setB - Set the B matrix of the PhasorSS object and update the isReal property
%   setC - Set the C matrix of the PhasorSS object and update the isReal property
%   setD - Set the D matrix of the PhasorSS object and update the isReal property
%   addOutput - Add additional output to the current output of the PhasorSS object
%   addInput - Add additional input to the current input of the PhasorSS object
%   expandBase - Expand the base of the PhasorSS object by a factor m
%   trunc - Truncate the PhasorSS object to h harmonics
%   neglect - Neglect the PhasorSS object below a threshold
%
%   % Manipulation methods
%   mtimes - Multiply two PhasorSS objects
%   feedback - Feedback connection of two PhasorSS objects
%   toLPVss - Convert the PhasorSS object to an LPV state-space system
%   toLTVss - Convert the PhasorSS object to an LTV state-space system
%   toSS - Convert the PhasorSS object to the appropriate ss object (wrapper for toLPVss and toLTVss)
%   evalAngle - Evaluate the PhasorSS object at a given angle
%
%   % Analysis methods
%   HmqBode - Compute the Bode plot of the PhasorSS object in the harmonic domain
%   hmqDcGain - Compute the DC gain of the PhasorSS object in the harmonic domain
%   toeplitzSS - Compute the Toeplitz Matrix State-Space object of the PhasorSS object
%
%   % Realness check methods
%   isreal - Check if the PhasorSS object is real
%   checkIfReal - Check if the PhasorArray in the PhasorSS object are real with a given tolerance, making the PhasorSS real
%   tolReal - Perform binary search to find the lowest tolerance that makes the system real
%   realify - Convert the PhasorSS object to a real-valued system
%
%   % LPV methods
%   setLPV - Set the LPV parameter of the PhasorSS object
%   checkp - Check the validity of the function p
%   setLTV - Set the PhasorSS object as a LTV
%
%   % Plot methods
%   plot - Plot the PhasorSS object
%   stem - Plot the stem of the PhasorSS object
%
%   % Utility methods
%   empty - Create an empty PhasorSS object
%   formatInputRange - Format the input range for inputHmRange
%   mergeStruct - Merge two structures
%
%   % Unsupported methods
%   cat - Concatenation is not supported for PhasorSS objects
%
% Example:
%   % Create a PhasorSS object with specified matrices
%   A = PhasorArray(cat(3,[0 1; -2 -3],[1/2 0; 0 0]),"isreal",true);
%   B = [0; 1];
%   C = [1 0];
%   D = 0;
%   T = 0.1;
%   Pss = PhasorSS(A, B, C, D, T, 'isReal', true); 
%   % Pss is a periodic real-valued state-space system with period 0.1 (A is a PhasorArray that define a cos for its (1,1) element)
%
%   % Set the LPV parameter for the phase
%   p = @(t,x,u) x(2); % LPV parameter for the phase
%   Pss = Pss.setLPV(p); 
%   % Pss is a "periodic" real-valued LPV state-space system where 
%   % periodic matrix are computed according to the angle defined by the LPV parameter p,
%   % here the phase is the second state, T is not used anymore.
%
%   % Perform a step input simulation
%   t = 0:0.01:10;
%   [y, tOut, x] = Pss.step(t);
%   % Plot the output response
%   plot(tOut, y);
%
% See also: PhasorArray


classdef PhasorSS < matlab.mixin.indexing.RedefinesParen & matlab.mixin.CustomDisplay

    
    properties (SetAccess = protected )
        A % State matrix (PhasorArray)
        B % Input matrix (PhasorArray)
        C % Output matrix (PhasorArray)
        D % Feedthrough matrix (PhasorArray)
        p % LPV parameter for phase used for periodicity
        isLPV % Logical flag indicating if the system is LPV or LTV
    end
    properties (SetAccess = public )
        T % Period of the system if set to LTV
        isReal % Logical flag indicating if the system is real
        StateName % Cell array of state names
        StateUnit % Cell array of state units
        InputName % Cell array of input names
        InputUnit % Cell array of input units
        InputGroup % Structure of input groups
        OutputName % Cell array of output names
        OutputUnit % Cell array of output units
        OutputGroup % Structure of output groups
        Name % Name of the system
        Notes % Cell array of notes
        UserData % User-defined data
    end
    properties (Hidden)
        providedD % hidden property to store the providedD and manage the empty D matrix case
        providedC % hidden property to store the providedC and manage the empty C matrix case (default is identity matrix)
    end

    methods
        function obj = PhasorSS(A,B,C,D,T,varg)
            %PHASORSS Construct an instance of this PhasorSS
            %   obj = PHASORSS(A,B,C,D,T,varg) creates an instance of the 
            %   PhasorSS class with the specified system matrices A, B, C, D, 
            %   and period T. Additional parameters can be specified 
            %   using the varg structure.
            %
            %   Inputs:
            %       A - State matrix
            %       B - Input matrix
            %       C - Output matrix (optional, default is identity matrix)
            %       D - Feedthrough matrix (optional, default is zero matrix)
            %       T - period of the system if set to LTV (optional, default is 1), if p is set and the system is LPV, T is not used as the phase is computed from p
            %       varg - Structure containing additional parameters:
            %           isReal - Logical flag indicating if the system is real (default is false)
            %           StateName - Cell array of state names (default is empty)
            %           StateUnit - Cell array of state units (default is empty)
            %           InputName - Cell array of input names (default is empty)
            %           InputUnit - Cell array of input units (default is empty)
            %           InputGroup - Structure of input groups (default is empty struct)
            %           OutputName - Cell array of output names (default is empty)
            %           OutputUnit - Cell array of output units (default is empty)
            %           OutputGroup - Structure of output groups (default is empty struct)
            %           Name - Name of the system (default is empty string)
            %           Notes - Cell array of notes (default is empty)
            %           UserData - User-defined data (default is empty struct)
            %           p - LPV parameter for phase used for periodicity (default is empty)
            %
            % Example:
            %   % Create a PhasorSS object with specified matrices
            %   A = PhasorArray(cat(3,[0 1; -2 -3],[1/2 0; 0 0]),"isreal",true);
            %   B = [0; 1];
            %   C = [1 0];
            %   D = 0;
            %   T = 0.1;
            %   Pss = PhasorSS(A, B, C, D, T, 'isReal', true); 
            %   % Pss is a periodic real-valued state-space system with period 0.1 (A is a PhasorArray that define a cos for its (1,1) element)
            %
            %   % Set the LPV parameter for the phase
            %   p = @(t,x,u) x(2); % LPV parameter for the phase
            %   Pss = Pss.setLPV(p); 
            %   % Pss is a "periodic" real-valued LPV state-space system where 
            %   % periodic matrix are computed according to the angle defined by the LPV parameter p,
            %   % here the phase is the second state, T is not used anymore.
            %
            %   % Perform a step input simulation
            %   t = 0:0.01:10;
            %   [y, tOut, x] = Pss.step(t);
            %   % Plot the output response
            %   plot(tOut, y);
            %
            % See also: PhasorArray
            
            arguments
                A
                B = []
                C = []
                D = []
                T = []
                varg.isReal (1,1) logical = false
                varg.StateName (1,:)  = {}
                varg.StateUnit (1,:)  = {}
                varg.InputName (1,:)  = {}
                varg.InputUnit (1,:)  = {}
                varg.InputGroup (1,:) struct = struct()
                varg.OutputName (1,:)  = {}
                varg.OutputUnit (1,:)  = {}
                varg.OutputGroup (1,:) struct = struct()
                varg.Name  = ''
                varg.Notes (1,:)  = {}
                varg.UserData = [];
                varg.p = []
                varg.verbose (1,1) logical = false
            end

            if nargin == 1
                if isa(A,'PhasorSS')
                    obj = A;
                    return
                end
                if isa(A,'ss') || isa(A,'tf')
                    obj = PhasorSS.fromSS(A);
                    return
                end
                D = A;
                A = [];
            end

            if D == 0 ; D=[]; end %if D is 0, set it to empty matrix to avoid confusion
            obj.providedC = C;
            obj.providedD = D;
            if isempty(C)
                C=eye(size(A,1));
                varg.OutputName = varg.StateName;
                varg.OutputUnit = varg.StateUnit      ;
            end
            if isempty(D) || (isscalar(D) && all(D ==0))
                D=zeros(size(C,1),size(B,2));
            end
            
            %if only D is provided, A is empty matrix, B is 0xsize(D,2), C is size(D,1)x0, T is 1
            if isempty(A)
                if isempty(B) && isempty(C)
                A = PhasorArray.empty(0);
                B = PhasorArray.empty(0,size(D,2));
                C = PhasorArray.empty(size(D,1),0);
                else
                    error('If A is specified empty, B and C must be empty')
                end
            end

            %finally, check that all dimension are consistent
            if size(A,1) ~= size(A,2)
                error('A must be a square matrix')
            end
            if size(A,1) ~= size(B,1)
                error('A and B must have the same number of rows')
            end
            if size(C,2) ~= size(A,1)
                error('C must have the same number of columns as A')
            end
            if size(D,1) ~= size(C,1)
                error('D must have the same number of rows as C')
            end
            if size(D,2) ~= size(B,2)
                error('D must have the same number of columns as B')
            end

            % fill input, output, state with u, y, x if empty
            if isempty(varg.InputName)
                varg.InputName = cellstr(strcat('u',num2str((1:size(B,2))')))';
            end
            if isempty(varg.OutputName)
                varg.OutputName = cellstr(strcat('y',num2str((1:size(C,1))')))';
            end
            if isempty(varg.StateName) && ~isempty(A)
                varg.StateName = cellstr(strcat('x',num2str((1:size(A,1))')))';
            end
            if isempty(varg.InputUnit)
                varg.InputUnit = repmat({''},1,size(B,2));
            end
            if isempty(varg.OutputUnit)
                varg.OutputUnit = repmat({''},1,size(C,1));
            end
            if isempty(varg.StateUnit) && ~isempty(A)
                varg.StateUnit = repmat({''},1,size(A,1));
            end

            %check that all dimensions are consistent
            if size(varg.InputName,2) ~= size(B,2)
                error('InputName must have the same number (%d) of rows as B (%d) ',size(varg.InputName,2),size(B,2))
            end
            if size(varg.InputUnit,2) ~= size(B,2)
                error('InputUnit must have the same number (%d) of rows as B (%d) ',size(varg.InputUnit,2),size(B,2))
            end
            if size(varg.OutputName,2) ~= size(C,1)
                error('OutputName must have the same number (%d) of rows as C (%d) ',size(varg.OutputName,2),size(C,1))
            end
            if size(varg.OutputUnit,2) ~= size(C,1)
                error('OutputUnit must have the same number (%d) of rows as C (%d) ',size(varg.OutputUnit,2),size(C,1))
            end
            if size(varg.StateName,2) ~= size(A,1)
                error('StateName must have the same number (%d) of rows as A (%d) ',size(varg.StateName,2),size(A,1))
            end
            if size(varg.StateUnit,2) ~= size(A,1)
                error('StateUnit must have the same number (%d) of rows as A (%d) ',size(varg.StateUnit,2),size(A,1))
            end


            obj.A = PhasorArray(A);
            obj.B = PhasorArray(B);
            obj.C = PhasorArray(C);
            obj.D = PhasorArray(D);
            obj.T = T;
            obj.isReal = varg.isReal;
            obj.StateName = cellstr(varg.StateName);
            obj.StateUnit = cellstr(varg.StateUnit);
            obj.InputName = cellstr(varg.InputName);
            obj.InputUnit = cellstr(varg.InputUnit);
            obj.InputGroup = varg.InputGroup;
            obj.OutputName = cellstr(varg.OutputName);
            obj.OutputUnit = cellstr(varg.OutputUnit);
            obj.OutputGroup = varg.OutputGroup;
            obj.Name = varg.Name;
            obj.Notes = cellstr(varg.Notes);
            obj.UserData = varg.UserData;

            obj = obj.checkIfReal(1e-12,varg.verbose);

            if ~isempty(varg.p)	
                obj = setLPV(obj,varg.p);
            else
                if isempty(obj.T)
                    obj.T = 1;
                    warning('phasorSS:noT','Using PHASORSS Constructor : Period T is empty, system is set as LTV with period 1')
                end
                obj = setLTV(obj);
            end
        end

        function obj = setLPV(obj,p)
            %SETLPV Set the LPV parameter of the PhasorSS object
            %   obj = SETLPV(obj,p) sets the LPV parameter of the PhasorSS object.
            %   If the parameter is empty, the system is set as an LTV system.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       p - LPV parameter (can be a function or indices)
            %
            %   Outputs:
            %       obj - Updated instance of the PhasorSS class
            %
            %   Example:
            %       obj = obj.setLPV(@(t,x,u) u(1));
            %       set the LPV phase parameter as the first input
            %
            %       obj = obj.setLPV(1);
            %       set the LPV phase parameter as the first state     
            %
            %   See also: SETLTV

            arguments
                obj
                p = []
            end
            if isempty(p)
                warning('phasorSS:notLPV','LPV parameter is empty, set as LTV system instead')
                obj = obj.setLTV();
                return
            end
            if isnumeric(p)
                ind = p;
                p = @(t,x,u) x(ind);
            end
            obj.p = p;
            checkp(obj);
            obj.isLPV = true;
        end
        function checkp(obj)
            % CHECKP Check the validity of the function p
            %   This method verifies that the function p is a 3-argument function
            %   compatible with the format [scalar, vector, vector] -> scalar.
            %   It ensures that p can accept a scalar and two vectors as inputs and
            %   returns a scalar output.
            %
            %   Usage:
            %       obj.checkp()
            %
            %   Example:
            %       obj.checkp();
            %
            %   See also: setLPV
            try
                ptest = obj.p(0,zeros(size(obj.A,1),1),zeros(size(obj.B,2),1));
                if ~isscalar(ptest)
                    error('LPV parameter p must return a scalar')
                end
            catch
                error('LPV parameter p must be a function @(t,x,u) returning a scalar')
            end
        end

        function obj = setLTV(obj,T)
            %SETLTV Set the PhasorSS object as a LTV 
            %   obj = SETLTV(obj,T) sets the PhasorSS object as a Linear Time-Varying (LTV) system, 
            %       optionally set the period T of obj.
            %   If the period T is empty in obj, a warning is displayed.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       T - Period of the system (default is empty)
            %
            %   Outputs:
            %       obj - Updated instance of the PhasorSS class
            %
            %   Example:
            %       T = 0.1;
            %       obj = obj.setLTV(T);
            %       Set the PhasorSS object as a LTV system with period T
            %
            %   See also: SETLPV
            arguments
                obj
                T  = []
            end
            if ~isempty(T)
                obj.T = T;
            end
            obj.isLPV = false;
            obj.p = [];
            if isempty(obj.T)
                warning('phasorSS:noT','Using SETLTV : Period T is empty, set a value to use the system as LTV')
            end
        end

        function out = isreal(obj)
            %ISREAL Check if the PhasorSS object is real
            %   out = ISREAL(obj) checks if the PhasorSS object is real.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       out - Boolean indicating if the system is real
            %
            %   See also: CHECKIFREAL

            obj = obj.checkIfReal();
            out = obj.isReal;
        end

        function obj = checkIfReal(obj,tol,verbose)
            %CHECKIFREAL Check if the PhasorArray in the PhasorSS object are real with a given tolerance, making the PhasorSS real 
            %   obj = CHECKIFREAL(obj,tol,verbose) checks if the PhasorArray in the PhasorSS object are real with a given tolerance.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       tol - Tolerance for checking if the system is real (default is 1e-12)
            %       verbose - Boolean indicating if warnings should be displayed (default is false)
            %
            %   Outputs:
            %       obj - Updated instance of the PhasorSS class with isReal property
            %
            %   See also: TOLREAL, ISREAL, REALIFY

            arguments
                obj
                tol (1,1) double = 1e-12
                verbose logical = false
            end
            [aReal] = isreal(obj.A,tol);
            [bReal] = isreal(obj.B,tol);
            [cReal] = isreal(obj.C,tol);
            [dReal] = isreal(obj.D,tol);

            if aReal && bReal && cReal && dReal
                obj.isReal = true;
                if ~obj.isReal || verbose
                    warning('phasorSS:appearsToBeReal','The PhasorSS object is real with a tolerance of %e. \n Manually set isReal to false to force complex valued computation (longer simulation time)',tol)
                end
                obj.isReal = true;
            else
                obj.isReal = false;
            end
        end

        function tol = tolReal(obj,tolstart,tolTol)
            %TOLREAL Perform binary search to find the lowest tolerance that makes the system real
            %   tol = TOLREAL(obj,tolstart,tolTol) performs a binary search to find the tolerance that makes the system real.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       tolstart - Starting tolerance for the binary search (default is 1e-4)
            %       tolTol - Tolerance for the binary search (default is 1e-20)
            %
            %   Outputs:
            %       tol - lowest tolerance that makes the system 
            %
            %   See also: CHECKIFREAL, ISREAL, REALIFY
            arguments
                obj
                tolstart (1,1) double = 1e-4
                tolTol (1,1) double = 1e-20
            end
            
            low = 0;
            high = tolstart;
            while high - low > tolTol
                mid = (high + low) / 2;
                obj = obj.checkIfReal(mid);
                if obj.isReal
                    high = mid;
                else
                    low = mid;
                end
            end
            tol = mid;
        end



        function objo = realify(obj)
            %REALIFY Convert the PhasorSS object to a real-valued system
            %   objo = REALIFY(obj) converts the PhasorSS object to a real-valued system.
            %
            %   See also: CHECKIFREAL, ISREAL, TOLREAL
            objo = obj;
            objo.setA(mreal(obj.A));
            objo.setB(mreal(obj.B));
            objo.setC(mreal(obj.C));
            objo.setD(mreal(obj.D));
            objo = objo.checkIfReal();
        end

        function obj = setA(obj,A,realTol)
            %SETA Set the A matrix of the PhasorSS object and update the isReal property if needed
            %   obj = SETA(obj,A,realTol) sets the A matrix of the PhasorSS object and updates the isReal property.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       A - State matrix (PhasorArray)
            %       realTol - Tolerance for checking if the system is real (default is 1e-12)
            %
            %   Outputs:
            %       obj - Updated instance of the PhasorSS class
            %
            %   See also: SETB, SETC, SETD, SETLPV, CHECKIFREAL
            arguments
                obj
                A
                realTol (1,1) double = 1e-12
            end
            %SETA Set the A matrix of the PhasorSS object and update the isReal property
            obj.A = PhasorArray(A);
            obj = obj.checkIfReal(realTol);


            if isempty(obj.providedC)
                obj = obj.setC([]);
            end
        end
        function obj = setB(obj,B,realTol,optargin)
            %SETB Set the B matrix of the PhasorSS object and update the isReal property
            %   obj = SETB(obj,B,realTol) sets the B matrix of the PhasorSS object and updates the isReal property.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       B - Input matrix (PhasorArray)
            %       realTol - Tolerance for checking if the system is real (default is 1e-12)
            %
            %   Outputs:
            %       obj - Updated instance of the PhasorSS class
            %
            %   See also: SETA, SETC, SETD, SETLPV, CHECKIFREAL
            arguments
                obj
                B
                realTol (1,1) double = 1e-12
                optargin.InputName (1,:) cell = {}
                optargin.InputUnit (1,:) cell = {}
                optargin.InputGroup (1,:) struct = struct()
            end
            %SETB Set the B matrix of the PhasorSS object and update the isReal property
            obj.B = PhasorArray(B);
            obj = obj.checkIfReal(realTol);

            %auto update D if not provided in case number of input change
            if isempty(obj.providedD)
                obj = obj.setD([]);
            end

            %update name and unit
            if isempty(optargin.InputName)
                obj.InputName = cellstr(strcat('u',num2str((1:size(B,2))')))';
            else
                obj.InputName = cellstr(optargin.InputName);
            end
            if isempty(optargin.InputUnit)
                obj.InputUnit = repmat({''},1,size(B,2));
            else
                obj.InputUnit = cellstr(optargin.InputUnit);
            end
            obj.InputGroup = optargin.InputGroup;
        end
        function obj = setC(obj,C,realTol,optargin)
            %SETC Set the C matrix of the PhasorSS object and update the isReal property
            %   obj = SETC(obj,C,realTol) sets the C matrix of the PhasorSS object and updates the isReal property.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       C - Output matrix (PhasorArray)
            %       realTol - Tolerance for checking if the system is real (default is 1e-12)
            %
            %   Outputs:
            %       obj - Updated instance of the PhasorSS class
            %
            %   See also: SETA, SETB, SETD, SETLPV, CHECKIFREAL
            arguments
                obj
                C
                realTol (1,1) double = 1e-12
                optargin.OutputName (1,:) cell = {}
                optargin.OutputUnit (1,:) cell = {}
                optargin.OutputGroup (1,:) struct = struct()
            end

            if isempty(C)
                C=eye(size(obj.A,1));
                obj.OutputName = obj.StateName;
                obj.OutputUnit = obj.StateUnit;
                obj.providedC = [];
            else
                obj.providedC = C;
            end
            obj.C = PhasorArray(C);
            if isempty(obj.providedD)
                obj.D=PhasorArray(zeros(size(obj.C,1),size(obj.B,2)));
            end
            obj = obj.checkIfReal(realTol);

            %update name and unit
            if isempty(optargin.OutputName)
                obj.OutputName = cellstr(strcat('y',num2str((1:size(C,1))')))';
            else
                obj.OutputName = cellstr(optargin.OutputName);
            end
            if isempty(optargin.OutputUnit)
                obj.OutputUnit = repmat({''},1,size(C,1));
            else
                obj.OutputUnit = cellstr(optargin.OutputUnit);
            end
            obj.OutputGroup = optargin.OutputGroup;
        end
        function obj = setD(obj,D,realTol)
            %SETD Set the D matrix of the PhasorSS object and update the isReal property
            %   obj = SETD(obj,D,realTol) sets the D matrix of the PhasorSS object and updates the isReal property.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       D - Feedthrough matrix (PhasorArray)
            %       realTol - Tolerance for checking if the system is real (default is 1e-12)
            %
            %   Outputs:
            %       obj - Updated instance of the PhasorSS class
            %
            %   See also: SETA, SETB, SETC, SETLPV, CHECKIFREAL
            arguments
                obj
                D
                realTol (1,1) double = 1e-12
            end

            if isempty(D)
                obj.D=PhasorArray(zeros(size(obj.C,1),size(obj.B,2)));
                obj.providedD = [];
                obj = obj.checkIfReal(realTol);
                return
            end

            obj.D = PhasorArray(D);
            obj = obj.checkIfReal(realTol);
            obj.providedD = D;
        end

        function nx = nx(obj)
            %NX Number of states of the PhasorSS object
            %   nx = NX(obj) returns the number of states of the PhasorSS object.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       nx - Number of states
            %
            %   See also: nU, nY
            nx = size(obj.A,1);
        end

        function nu = nu(obj)
            %NU Number of inputs of the PhasorSS object
            %   nu = NU(obj) returns the number of inputs of the PhasorSS object.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       nu - Number of inputs
            %
            %   See also: nX, nY
            nu = size(obj.B,2);
        end

        function ny = ny(obj)
            %NY Number of outputs of the PhasorSS object
            %   ny = NY(obj) returns the number of outputs of the PhasorSS object.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       ny - Number of outputs
            %
            %   See also: nX, nU
            ny = size(obj.C,1);
        end

        function obj = addOutput(obj,newC,newD,varg)
            %ADDOUTPUT Add additionnal output to the current output of the PhasorSS object
            %  obj = ADDOUTPUT(obj,newC,newD,varg) adds additionnal output to the current output of the PhasorSS object.
            %  If newC is empty, newD must be provided and vice versa.
            %  If newC and newD are empty, an error is thrown.
            %
            %  Inputs:
            %      obj - Instance of the PhasorSS class
            %      newC - Output matrix (PhasorArray)
            %      newD - Feedthrough matrix (PhasorArray)
            %      varg - Structure containing additional parameters:
            %          OutputName - Cell array of output names (default is empty)
            %          OutputUnit - Cell array of output units (default is empty)
            %          OutputGroup - Structure of output groups (default is empty struct)
            %
            %   see also: ADDINPUT, SETC, SETD
            arguments
                obj
                newC = []
                newD = []
                varg.OutputName (1,:) cell = {}
                varg.OutputUnit (1,:) cell = {}
                varg.OutputGroup (1,:) struct = struct()
            end
            if isempty(newC)
                if isempty(newD)
                    error('newC and newD cannot be both empty')
                end
                newC = zeros(size(newD,1),size(obj.A,1));
            end
            if isempty(newD)
                newD = zeros(size(newC,1),size(obj.B,2));
            end

            if size(newC,2) ~= size(obj.A,1)
                error('newC must have the same number of columns as A')
            end
            if size(newD,1) ~= size(newC,1)
                error('newD must have the same number of rows as newC')
            end
            if size(newD,2) ~= size(obj.B,2)
                error('newD must have the same number of columns as B')
            end

            obj.C = [obj.C;PhasorArray(newC)];
            obj.D = [obj.D;PhasorArray(newD)];

            %update name and unit
            if isempty(varg.OutputName)
                obj.OutputName = [obj.OutputName cellstr(strcat('y',num2str((size(obj.C,1)-size(newC,1)+1:size(obj.C,1))')))];
            else
                obj.OutputName = [obj.OutputName cellstr(varg.OutputName)];
            end
            if isempty(varg.OutputUnit)
                obj.OutputUnit = [obj.OutputUnit repmat({''},1,size(newC,1))];
            else
                obj.OutputUnit = [obj.OutputUnit cellstr(varg.OutputUnit)];
            end

            obj.OutputGroup = mergeStruct(obj.OutputGroup,varg.OutputGroup);
            
            obj.providedC = obj.C;
            obj.providedD = obj.D;

        end

        function obj = addInput(obj,newB,newD,varg)
            %ADDINPUT Add additionnal input to the current input of the PhasorSS object
            %  obj = ADDINPUT(obj,newB,newD,varg) adds additionnal input to the current input of the PhasorSS object.
            %  If newB is empty, newD must be provided and vice versa.
            %  If newB and newD are empty, an error is thrown.
            %
            %  Inputs:
            %      obj - Instance of the PhasorSS class
            %      newB - Input matrix (PhasorArray)
            %      newD - Feedthrough matrix (PhasorArray)
            %      varg - Structure containing additional parameters:
            %          InputName - Cell array of input names (default is empty)
            %          InputUnit - Cell array of input units (default is empty)
            %          InputGroup - Structure of input groups (default is empty struct)
            %
            %  Outputs:
            %      obj - Updated instance of the PhasorSS class
            %
            %  See also: ADDOUTPUT, SETB, SETD
            arguments
                obj
                newB = []
                newD = []
                varg.InputName (1,:) cell = {}
                varg.InputUnit (1,:) cell = {}
                varg.InputGroup (1,:) struct = struct()
            end
            if isempty(newB)
                if isempty(newD)
                    error('newB and newD cannot be both empty')
                end
                newB = zeros(size(obj.A,1),size(newD,2));
            end
            if isempty(newD)
                newD = zeros(size(obj.C,1),size(newB,2));
            end

            if size(newB,1) ~= size(obj.A,1)
                error('newB must have the same number of rows as A')
            end
            if size(newD,2) ~= size(newB,2)
                error('newD must have the same number of columns as newB')
            end
            if size(newD,1) ~= size(obj.C,1)
                error('newD must have the same number of rows as C')
            end

            obj.B = [obj.B PhasorArray(newB)];
            obj.D = [obj.D PhasorArray(newD)];

            %update name and unit
            if isempty(varg.InputName)
                obj.InputName = [obj.InputName cellstr(strcat('u',num2str((size(obj.B,2)-size(newB,2)+1:size(obj.B,2))')))];
            else
                obj.InputName = [obj.InputName cellstr(varg.InputName)];
            end
            if isempty(varg.InputUnit)
                obj.InputUnit = [obj.InputUnit repmat({''},1,size(newB,2))];
            else
                obj.InputUnit = [obj.InputUnit cellstr(varg.InputUnit)];
            end
            obj.InputGroup = mergeStruct(obj.InputGroup,varg.InputGroup);
        end

        function sys = toLPVss(obj)
            %TOLPVSS Convert the PhasorSS object to an LPV state-space system
            %   sys = TOLPVSS(obj) converts the PhasorSS object to a Linear Parameter-Varying (LPV) state-space system.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       sys - LPV state-space system (matlab native object)
            %   Example:
            %       sys = obj.toLPVss();
            %       Convert the PhasorSS object to an LPV state-space system
            %
            %   See also: TOLTVSS, TOSS
            if obj.isReal
                disp('Real valued LPP ss')
                Acs = obj.A.SinCosForm();
                Bcs = obj.B.SinCosForm();
                Ccs = obj.C.SinCosForm();
                Dcs = obj.D.SinCosForm();

                ltvFun = @(t,p) lpvFunSC(t,p,Acs,Bcs,Ccs,Dcs);

            else
                disp('Complex valued LPP ss')
                A = obj.A;
                B = obj.B;
                C = obj.C;
                D = obj.D;

                ltvFun = @(t,p) lpvFunP(t,p,A,B,C,D);

            end
            sys = lpvss("phase",ltvFun, ...
                'StateName',obj.StateName,'StateUnit',obj.StateUnit, ...
                'InputName',obj.InputName,'InputUnit',obj.InputUnit,'InputGroup',obj.InputGroup, ...
                'OutputName',obj.OutputName,'OutputUnit',obj.OutputUnit,'OutputGroup',obj.OutputGroup, ...
                'Name',obj.Name,'Notes',obj.Notes,'UserData',obj.UserData);


            function [A,B,C,D,E,dx0,x0,u0,y0,Delays] = lpvFunSC(~,p,Acs,Bcs,Ccs,Dcs)
                %LPVFUNSC LPV function for the PhasorSS object in the SinCosForm
                E = [];
                dx0 = [];
                x0 = [];
                u0 = [];
                y0 = [];
                Delays = [];

                A=sincos2time(Acs,p);
                B=sincos2time(Bcs,p);
                C=sincos2time(Ccs,p);
                D=sincos2time(Dcs,p);
            end
            function [A,B,C,D,E,dx0,x0,u0,y0,Delays] = lpvFunP(~,p,A,B,C,D)
                %LPVFUNP LPV function for the PhasorSS object in the PhasorForm
                E = [];
                dx0 = [];
                x0 = [];
                u0 = [];
                y0 = [];
                Delays = [];

                A=phasor2time(A,p);
                B=phasor2time(B,p);
                C=phasor2time(C,p);
                D=phasor2time(D,p);
            end


        end


        function sys = toLTVss(obj)
            %TOLTVSS Convert the PhasorSS object to an LTV state-space system
            %   sys = TOLTVSS(obj) converts the PhasorSS object to a Linear Time-Varying (LTV) state-space system.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       sys - LTV state-space system (matlab native object)
            %   Example:
            %       sys = obj.toLTVss();
            %       Convert the PhasorSS object to an LTV state-space system
            %
            %   See also: TOLPVSS, TOSS
            if isempty(obj.T)
                error('phasorSS:noT','Using TOLTVSS : Period T is empty, set a value to use the system as LTV')
            end

            if obj.isReal
                disp('Real valued LTP ss')
                Acs = obj.A.SinCosForm();
                Bcs = obj.B.SinCosForm();
                Ccs = obj.C.SinCosForm();
                Dcs = obj.D.SinCosForm();

                ltvFun = @(t) ltvFunSC(t,Acs,Bcs,Ccs,Dcs);

            else
                disp('Complex valued LTP ss')
                ltvFun = @(t) ltvFunP(t);
            end
            sys = ltvss(ltvFun,'StateName',obj.StateName,'StateUnit',obj.StateUnit,'InputName',obj.InputName,'InputUnit',obj.InputUnit,'OutputName',obj.OutputName,'OutputUnit',obj.OutputUnit,'InputGroup',obj.InputGroup,'OutputGroup',obj.OutputGroup,'Name',obj.Name,'Notes',obj.Notes,'UserData',obj.UserData);


            function [A,B,C,D,E,dx0,x0,u0,y0,Delay] = ltvFunSC(t,Acs,Bcs,Ccs,Dcs)
                %LTVFUNSC LTV function for the PhasorSS object in the SinCosForm

                E      = [];
                dx0    = [];
                x0     = [];
                u0     = [];
                y0     = [];
                Delay  = [];
                
                %Delay.Output = NaN(1,size(Ccs,1))';
                %Delay.Input  = NaN(1,size(Bcs,2))';

                A = double(sincos2time(Acs,t/obj.T*2*pi));
                B = double(sincos2time(Bcs,t/obj.T*2*pi));
                C = double(sincos2time(Ccs,t/obj.T*2*pi));
                D = double(sincos2time(Dcs,t/obj.T*2*pi));
            end
            function [A,B,C,D,E,dx0,x0,u0,y0,Delays] = ltvFunP(t)
                %LTVFUNP LTV function for the PhasorSS object in the PhasorForm

                E = [];
                dx0 = [];
                x0 = [];
                u0 = [];
                y0 = [];
                Delays = [];

                A=double(phasor2time(obj.A.value,t/obj.T*2*pi));
                B=double(phasor2time(obj.B.value,t/obj.T*2*pi));
                C=double(phasor2time(obj.C.value,t/obj.T*2*pi));
                D=double(phasor2time(obj.D.value,t/obj.T*2*pi));
            end


        end

        function sys = toSS(obj)
            %TOSS Convert the PhasorSS object to the appropriate ss object (wrapper for TOLPVSS and TOLTVSS)
            %   sys = TOSS(obj) converts the PhasorSS object to the appropriate state-space object in native matlab format.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       sys - State-space object (matlab native object)
            %   Example:
            %       sys = obj.toSS();
            %       Convert the PhasorSS object to the appropriate state-space object
            %
            %   See also: TOLPVSS, TOLTVSS
            if obj.isLPV
                sys = obj.toLPVss();
            else
                sys = obj.toLTVss();
            end
        end

        function varargout = initial(sys,x0,t)
            %INITIAL Perform initial condition simulation of the PhasorSS object
            %   INITIAL(sys,x0,t) performs initial condition simulation of the PhasorSS object.
            %
            %   Inputs:
            %       sys - Instance of the PhasorSS class
            %       x0 - Initial condition, default value ones(nx,1)
            %       t - Time vector, default value 1
            %
            %   Outputs (in order):
            %       y - Output response
            %       tOut - Output time vector
            %       x - State response
            %       pOut - Parameter response (if LPV system)
            %
            %   Example:
            %       initial(sys,x0,t);
            %       Perform initial condition simulation of the PhasorSS object
            %
            %   See also: STEP, IMPULSE, LSIM
            arguments
                sys
                x0 = ones(size(sys.A,1),1);
                t  = 1;
            end
            sysltv = sys.toSS();
            if sys.isLPV
                args = {sysltv,{x0,sys.p(0,x0,0)},t,sys.p};
            else
                args = {sysltv,x0,t};
            end
            [varargout{1:nargout}] = initial(args{:});            
        end

        function [y,tOut,x,pOut] = step(sys,t,stepOpt)
            %STEP Perform step input simulation of the PhasorSS object
            %   STEP(sys,t) performs step input simulation of the PhasorSS object.
            %   note that for LPV system, initial condition is set such that p(0,x0,0) = 0
            %
            %   Inputs:
            %       sys - Instance of the PhasorSS class
            %       t - Time vector
            %
            %   Example:
            %       step(sys,t);
            %       Perform step input simulation of the PhasorSS object
            %
            %   See also: INITIAL, IMPULSE, LSIM
            arguments
                sys
                t
                stepOpt = []
            end
            sysltv = sys.toSS();
            if sys.isLPV
                respOpt = RespConfig(InitialState = 'x0',InitialParameter=0);
                %merge stepOpt with respOpt
                if ~isempty(stepOpt)
                    fn = fieldnames(stepOpt);
                    for i = 1:numel(fn)
                        respOpt.(fn{i}) = stepOpt.(fn{i});
                    end                    
                end
                if nargout > 0
                    [y,tOut,x,~,pOut] = step(sysltv,t,sys.p,respOpt);
                else
                    step(sysltv,t,sys.p,respOpt);
                end
            else
                if nargout > 0
                [y,tOut,x] = step(sysltv,t,stepOpt);
                    if nargout > 3
                        pOut = tOut*2*pi/sys.T;
                        warning('phasorSS:noP','No parameter output for non LPV system, pOut is set to t*2*pi/T')
                    end
                else
                    step(sysltv,t,stepOpt);
                end
            end
        end

        function [yout,t,x] = impulse(sys,t)
            %IMPULSE Perform impulse input simulation of the PhasorSS object
            %   IMPULSE(sys,t) performs impulse input simulation of the PhasorSS object.
            %
            %   Inputs:
            %       sys - Instance of the PhasorSS class
            %       t - Time vector
            %
            %   Example:
            %       impulse(sys,t);
            %       Perform impulse input simulation of the PhasorSS object
            %
            %   See also: INITIAL, STEP, LSIM
            arguments
                sys
                t
            end
            sysltv = sys.toSS();
            if sys.isLPV
                respOpt = RespConfig( InitialState = 'x0',InitialParameter=0);

                if nargout > 0
                [yout,t,x] = impulse(sysltv,t,sys.p,respOpt);
                else
                impulse(sysltv,t,sys.p,respOpt)
                end
            else
                if nargout > 0
                [yout,t,x] = impulse(sysltv,t);
                else
                impulse(sysltv,t)
                end
            end
        end

        function [y,tOut,x,pOut] = lsim(sys,t,u,x0)
            %LSIM Perform simulation of the PhasorSS object, with input history (u,t) and initial condition x0
            %   [y,tOut,x,pOut] = LSIM(sys,t,u,x0) performs simulation of the PhasorSS object, with input history (u,t) and initial condition x0.
            %
            %   Inputs:
            %       sys - Instance of the PhasorSS class
            %       t - Time vector
            %       u - Input history
            %       x0 - Initial condition (default is zeros)
            %
            %   Outputs:
            %       y - Output response
            %       tOut - Output time vector
            %       x - State response
            %       pOut - Parameter response (if LPV system)
            %
            %   Example:
            %       [y,tOut,x,pOut] = lsim(sys,t,u,x0);
            %       Perform simulation of the PhasorSS object
            %
            %   See also: LSIMPLOT, STEP, IMPULSE

            arguments
                sys
                t
                u
                x0 = zeros(size(sys.A,1),1)
            end
            
            if isa(u,"PhasorArray")
                sys = sys.setB(sys.B*diag(u));
                sys = sys.setD(sys.D*diag(u));
                u = ones(size(u,1),numel(t));
            end

            sys = sys.checkIfReal(1e-12,true);
            sysltv = sys.toSS();
            if nargout == 0
                if sys.isLPV
                    lsim(sysltv,u,t,x0,sys.p)
                else
                    lsim(sysltv,u,t,x0)
                end
                return
            end
            if sys.isLPV
                [y,tOut,x,pOut] = lsim(sysltv,u,t,x0,sys.p);
            else
                [y,tOut,x] = lsim(sysltv,u,t,x0);
            end
        end

        function h = lsimplot(sys,t,u,x0,plotoptions)
            %LSIMPLOT Perform simulation of the PhasorSS object, with input history (u,t) and initial condition x0
            %   h = LSIMPLOT(sys,t,u,x0,plotoptions) performs simulation of the PhasorSS object, with input history (u,t) and initial condition x0.
            %
            %   Inputs:
            %       sys - Instance of the PhasorSS class
            %       t - Time vector
            %       u - Input history
            %       x0 - Initial condition (default is zeros)
            %       plotoptions - Plot options (default is empty)
            %           See lsimplot for more details on the plot options
            %
            %   Outputs:
            %       h - Plot handle
            %
            %   Example:
            %       t = 0:0.01:10;
            %       u = sin(t);
            %       h = lsimplot(sys,t,u,x0,plotoptions);
            %       Perform simulation of the PhasorSS object
            %
            %   See also: LSIM, STEP, IMPULSE
            arguments
                sys
                t
                u
                x0 = zeros(size(sys.A,1),1)
                plotoptions = []
            end
            sysltv = sys.toSS();

            if sys.isLPV
                h = lsimplot(sysltv,u,t,x0,sys.p,plotoptions);
            else
                h = lsimplot(sysltv,u,t,x0,plotoptions);
            end
        end


        function [y,tOut,x,pOut] = stepu(sys,t,u,stepOption)
            %STEPU Perform step input simulation of the PhasorSS object with a periodic input u
            %   [y, t, x, pOut] = STEPU(sys, t, u, stepOption) performs step input simulation of the PhasorSS object with a periodic input u.
            %   Note that for LPV system initial condition is computed so that p(0,x0,0) = 0.
            %
            %   Inputs:
            %       sys - Instance of the PhasorSS class
            %       t - Time vector
            %       u - Periodic input (PhasorArray)
            %       stepOption - (Optional) Structure containing additional step options compatible with RespConfig
            %
            %   Outputs:
            %       y - Output response
            %       t - Output time vector
            %       x - State response
            %       pOut - Parameter response (if LPV system)
            %
            %   Example:
            %       t = 0:0.01:10;
            %       u = PhasorArray(cat(3, ones(size(sys.B, 2), 1), zeros(size(sys.B, 2), 1), ones(size(sys.B, 2), 1))); % cosine input
            %       [y, t, x, pOut] = stepu(sys, t, u);
            %       % Perform step input simulation of the PhasorSS object with a periodic input u
            %
            %   See also: STEP, IMPULSE, LSIM, RespConfig
            arguments
                sys
                t
                u = eye(size(sys.B,2));
                stepOption = []
            end
            U = PhasorArray(u);
            BU = sys.B*diag(U);
            %avoid the warning check when the system is real
            warning off phasorSS:appearsToBeReal
            %set the B matrix of the system with the periodic input u
            sysBU = setB(sys,BU,"InputName",sys.InputName,"InputUnit",sys.InputUnit,"InputGroup",sys.InputGroup);
            %turn on the warning on phasorSS:appearsToBeReal
            warning on phasorSS:appearsToBeReal


            sysltv = sysBU.toSS();
            if sys.isLPV
                respOpt = RespConfig(InitialState = 'x0',InitialParameter=0);
                %merge stepOption with respOpt
                if ~isempty(stepOption)
                    fn = fieldnames(stepOption);
                    for i = 1:numel(fn)
                        respOpt.(fn{i}) = stepOption.(fn{i});
                    end                    
                end
                if nargout > 0
                    [y,tOut,x,~,pOut] = step(sysltv,t,sys.p,respOpt);
                else
                    step(sysltv,t,sys.p,respOpt)
                end
            else
                if nargout > 0
                    [y,tOut,x] = step(sysltv,t,stepOption);
                    if nargout > 3
                        pOut = tOut*2*pi/sys.T;
                        warning('phasorSS:noP','No parameter output for non LPV system, pOut is set to t*2*pi/T')
                    end
                else
                step(sysltv,t,stepOption)
                end
            end
        end
        function SSo = expandBase(SSp,m)
            %EXPANDBASE Expand the base of the PhasorSS object by a factor m
            %   SSo = EXPANDBASE(SSp, m) expands the base of the PhasorSS object SSp by a factor m.
            %
            %   Inputs:
            %       SSp - Instance of the PhasorSS class
            %       m - Factor to expand the base (integer)
            %
            %   Outputs:
            %       SSo - Expanded instance of the PhasorSS class
            %
            %   Example:
            %       SSo = expandBase(SSp, 2);
            %       % Expands the base of the PhasorSS object SSp by a factor of 2
            %
            %   See also: REDUCEBASE
            SSo = SSp;
            Ao = SSp.A.expandBase(m);
            SSo=SSo.setA(Ao);
            Bo = SSp.B.expandBase(m);
            SSo=SSo.setB(Bo);

            if ~isempty(SSp.providedC)
                Co = SSp.C.expandBase(m);
                SSo=SSo.setC(Co);
            else
                SSo=SSo.setC([]);
            end
            if~isempty(SSp.providedD)
                Do = SSp.D.expandBase(m);
                SSo=SSo.setD(Do);
            else
                SSo=SSo.setD([]);
            end

            SSo.T = SSo.T*m;

            SSo.p = @(t,x,u) SSo.p(t,x,u)/m;

            if SSo.isLPV
                warning('Base changed, please check the p function')
            end
        end

        function [Hsel,H, ssTB] = HmqBode(o1, h, T, h_compute, varg)
            %HMQBODE Compute the Bode plot of the PhasorSS object in the harmonic domain
            %   [mag, phase, freq, ssTB] = HMQBODE(o1, h, T, h_compute) computes the Bode plot of the PhasorSS object in the harmonic domain.
            %
            %   Inputs:
            %       o1 - Instance of the PhasorSS class
            %       h - Number of harmonics to consider (truncation) (integer)
            %       T - Period of the system (default is 1) to perform evaluation in the harmonic domain
            %       h_compute - Number of harmonics to compute the Toeplitz Matrix (default is h)
            %       Name-Value Pair Arguments:
            %           inputHmRange - Range of harmonics to consider for the input (default is {':',-h:h})
            %           outputHmRange - Range of harmonics to consider for the output (default is {':',-h:h})
            %           freqRange - Frequency range for the Bode plot (default is logspace(-2, 2, 100))
            %
            %   Specifying both h and h_compute allows to compute the toeplitz matrix and Bode plot with greater precision
            %   while plotting the Bode plot with only h harmonics. Note that h_compute must be greater or equal to h.
            %
            %   Outputs:
            %       mag - Magnitude of the Bode plot
            %       phase - Phase of the Bode plot
            %       freq - Frequency vector
            %       ssTB - Toeplitz Matrix State-Space object
            %
            %   Example:
            %       [mag, phase, freq, ssTB] = HmqBode(o1, 5, 1, 10, 'freqRange', logspace(-1, 2, 50));
            %
            %   See also: HMQDCGAIN, BODE
        
            arguments
                o1
                h
                T = 1
                h_compute = h
                varg.inputHmRange = {':', -h:h}
                varg.outputHmRange = {':', -h:h}
                varg.freqRange = logspace(-2, 2, 100)
            end
        
            %format input and output range
            outputHmRange = formatInputRange(o1, varg.outputHmRange, 'output');
            inputHmRange = formatInputRange(o1, varg.inputHmRange, 'input');

            % Compute the Toeplitz Matrix State-Space object
            ssTB = o1.toeplitzSS(h_compute,T);

            [outputHmRange,outIdx] = formatIndices(outputHmRange, size(ssTB.C, 1),h_compute);
            [inputHmRange,inIdx] = formatIndices(inputHmRange, size(ssTB.B, 2),h_compute);

        
            % Frequency range for the Bode plot
            freq = varg.freqRange;
            
            H = freqresp(ssTB, 2*pi*freq);

            Hsel = H(outIdx, inIdx, :);
        
            % Plot the Bode plot
            figure;
            for ii = 1:size(Hsel, 1)
                for jj = 1:size(Hsel, 2)
                    subplot(size(Hsel, 1), size(Hsel, 2), (ii-1)*size(Hsel, 2) + jj);
                    mag = abs(squeeze(Hsel(ii, jj, :)));
                    phase = angle(squeeze(Hsel(ii, jj, :))) * 180 / pi;
                    semilogx(freq, 20*log10(mag));
                    yyaxis right;
                    semilogx(freq, phase);
                    title(['Input ' ssTB.InputName{inIdx(jj)} ' to Output ' ssTB.OutputName{outIdx(ii)}]);
                    xlabel('Frequency (Hz)');
                    ylabel('Magnitude (dB) / Phase (degrees)');
                    grid on;
                end
            end
        end
        
        function [dcGain,ssTB] = hmqDcGain(o1,h,T,h_compute,varg)
            %HMQDCGAIN Compute or plot the DC gain of the PhasorSS object in the harmonic domaine
            %   [dcGain,ssTB] = HMQDCGAIN(o1,h,T) computes the DC gain of the PhasorSS object in the harmonic domaine.
            %   HMQDCGAIN(o1,h,T) plot the DC gain of the PhasorSS object using barsurf.
            %   Inputs:
            %       o1 - Instance of the PhasorSS class
            %       h - Number of harmonics to consider (truncation) (integer)
            %       T - Period of the system (default is 1) to perform evaluation in the harmonic domain
            %       h_compute - Number of harmonics to compute the Toeplitz Matrix (default is h)
            %       Name-Value Pair Arguments:
            %           inputHmRange - Range of harmonics to consider for the input (default is {':',-h:h})
            %           outputHmRange - Range of harmonics to consider for the output (default is {':',-h:h})
            %
            %   Specifying both h and h_compute allows to compute the toeplitz matrix and dc gain with greater precision
            %   while plotting the dc gain with only h harmonics. Note that h_compute must be greater or equal to h.
            %
            %   Outputs:
            %       dcGain - DC gain of the system
            %       ssTB - Toeplitz Matrix State-Space object
            %
            %   For greater possibilities in term of plotting, use the outputs dcGain and ssTB to extract the blocks and labels as in the example below: 
            %       [reducedDCGain,outputLabels,inputLabels] = extractBlocksAndLabels(dcGain,outputHmRange,inputHmRange,h_compute,ssTB.OutputName,ssTB.InputName);
            %       barsurf(reducedDCGain,'xticklabel',inputLabels,'yticklabel',outputLabels);
            %       Replace outputHmRange (default {':',-h:h}) by the desired range of harmonics to plot, for example
            %       {1,0,2,'-h:h'} to plot phasor 0 of first output and harmonics -h to h of second output
            %       The same for inputHmRange (default {':',-h:h})            %
            %
            %   Example:
            %       [dcGain,ssTB] = hmqDcGain(o1,h,T);
            %       Compute the DC gain of the PhasorSS object in the harmonic domaine
            %
            %   See also: TOEPLITZSS
            arguments
                o1
                h (1,1) double
                T double
                h_compute = []
                varg.inputHmRange = {':',-h:h};
                varg.outputHmRange = {':',-h:h};
            end
            outputHmRange = formatInputRange(o1, varg.outputHmRange,'output');
            inputHmRange  = formatInputRange(o1, varg.inputHmRange,'input');
            if isempty(h_compute)
                h_compute = h;
            end 
            assert(h_compute >= h, 'h_compute must be greater or equal to h')
            assert(isscalar(T), 'T must be a scalar')

            %compute Toeplitz Matrix of the system
            ssTB1 = o1.toeplitzSS(h_compute,T);
            %compute the DC gain of the system
            dcGain1 = dcgain(ssTB1);
            [reducedDCGain,outputLabels,inputLabels] = extractBlocksAndLabels(dcGain1,outputHmRange,inputHmRange,h_compute,ssTB1.OutputName,ssTB1.InputName);
            if nargout == 0
                %no output, display the DC gain as barsurf
                barsurf(reducedDCGain,'xticklabel',inputLabels,'yticklabel',outputLabels);
            else
                dcGain = reducedDCGain;
                ssTB = ssTB1;
            end
        end
        function formattedRange = formatInputRange(o1, inputRange,InOrOut)
            %FORMATINPUTRANGE Format the input range for inputHmRange
            %   formattedRange = FORMATINPUTRANGE(o1, inputRange) formats the input range for inputHmRange.
            %   If the input is a string, it matches it with an input name of o1 and replaces it with the index.
            %   If the input is already an integer, a vector, a double, or ':', it leaves it as it is.
            %   If the input is a cell containing strings, it matches every string with an index and replaces the cell with the list of indices.
            %
            %   Inputs:
            %       o1 - Instance of the PhasorSS class
            %       inputRange - Input range to format
            %
            %   Outputs:
            %       formattedRange - Formatted input range

            formattedRange = inputRange;
            if strcmp(InOrOut,'output')
                inputNames = o1.OutputName;
            else
                inputNames = o1.InputName;
            end

            for i = 1:2:length(inputRange)
                input = inputRange{i};
                if ischar(input) && ~strcmp(input,':')
                    idx = find(strcmp(inputNames, input));
                    if isempty(idx)
                        error('Input name "%s" not found in the PhasorSS object.', input);
                    end
                    formattedRange{i} = idx;
                elseif iscell(input)
                    idxList = [];
                    for j = 1:length(input)
                        idx = find(strcmp(inputNames, input{j}));
                        if isempty(idx)
                            error('Input name "%s" not found in the PhasorSS object.', input{j});
                        end
                        idxList = [idxList, idx];
                    end
                    formattedRange{i} = idxList;
                end
            end
        end


        function toeplitzSS = toeplitzSS(o1,h,T)
            %TOEPLITZSS Compute the Toeplitz Matrix State-Space object of the PhasorSS object
            %   toeplitzSS = TOEPLITZSS(o1,h,T) computes the Toeplitz Matrix State-Space object of the PhasorSS object.
            %
            %   Inputs:
            %       o1 - Instance of the PhasorSS class
            %       h - Number of harmonics to consider (truncation) (integer)
            %       T - Period of the system (default is 1) to perform evaluation in the harmonic domain
            %
            %   Outputs:
            %       toeplitzSS - Toeplitz Matrix State-Space object
            %
            %   Example:
            %       h = 5;
            %       T = 0.1;
            %       toeplitzSS = toeplitzSS(o1,h,T);
            %       Compute the Toeplitz Matrix State-Space object of the PhasorSS object
            %
            %   See also: HMQDCGAIN, LTVSS, LPVSS
            arguments
            o1
            h (1,1) double = max([o1.A.h,o1.B.h,o1.C.h,o1.D.h]);
            T  = []
            end
                %compute Toeplitz Matrix of the system
                A = o1.A;
                B = o1.B;
                C = o1.C;
                D = o1.D;

                ATB = A.T_tb(h);
                CTB = C.T_tb(h);
                DTB = D.T_tb(h);
                BTB = B.T_tb(h);

                hmqStateName={};
                hmqInputName={};
                hmqOutputName={};


                for ii = 1:numel(o1.StateName)
                    temp = "〈"+o1.StateName{ii}+"〉_{"+string(-h:h)+"}";
                    hmqStateName=[hmqStateName temp];
                end

                for ii = 1:numel(o1.InputName)
                    temp = "〈"+o1.InputName{ii}+"〉_{"+string(-h:h)+"}";
                    hmqInputName=[hmqInputName temp];
                end

                for ii = 1:numel(o1.OutputName)
                    temp = "〈"+o1.OutputName{ii}+"〉_{"+string(-h:h)+"}";
                    hmqOutputName=[hmqOutputName temp];
                end


                if ~isempty(T)
                    N = N_tb(A,h,T);

                    toeplitzSS = ss(ATB-N,BTB,CTB,DTB,'StateName',hmqStateName,'InputName',hmqInputName,'OutputName',hmqOutputName);
                else
                    toeplitzSS = @(T) ss(ATB-N_tb(A,h,T),BTB,CTB,DTB,'StateName',hmqStateName,'InputName',hmqInputName,'OutputName',hmqOutputName);
                end

        end

        function bool = isStatic(obj)
            %ISSTATIC Check if the PhasorSS object is static
            %   bool = ISSTATIC(obj) checks if the PhasorSS object is static.
            %  A static system is a system where the state matrix A is empty.

            %check if the system is static
            bool = isempty(obj.A);
        end

        function objOut = mtimes(obj,obj2)
            %MTIMES Multiply two PhasorSS objects
            %   MTIMES(obj,obj2) multiplies two PhasorSS objects.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       obj2 - Instance of the PhasorSS , or PhasorArray class, or matrix
            %

            %PhasorSS by PhasorSS product is still WIP :
                %- if phasorSS is only a D matrix (static input output, A,B,C are empty) and phasorSS2 is a full phasorSS, it should be possible
                %- the reverse should be possible too
                %- if both are full phasorSS, it should be possible to multiply them, but still and undefined yet

            % check dimension of D and the phasorSS to be sure they match
            % propagate correctly the input/output name and unit

            if ~isa(obj,'PhasorSS')
                if isscalar(obj)
                    obj = PhasorArray(obj*eye(size(obj2,1)));
                elseif isa(obj,'double')
                    obj = PhasorArray(obj);
                elseif isa(obj,'PhasorArray')
                else
                    error('The first argument must be a PhasorSS, a PhasorArray, a matrix or a scalar')
                end
                %convert it to a static PhasorSS
                obj = PhasorSS([],[],[],obj);
            end

            if ~isa(obj2,'PhasorSS')
                if isscalar(obj2)
                    obj2 = PhasorArray(obj2*eye(size(obj.D,2)));
                elseif isa(obj2,'double')
                    obj2 = PhasorArray(obj2);
                elseif isa(obj2,'PhasorArray')
                else
                    error('The second argument must be a PhasorSS, a PhasorArray, a matrix or a scalar')
                end
                %convert it to a static PhasorSS
                obj2 = PhasorSS([],[],[],obj2);
            end	
            % check if both are PhasorSS
            if isa(obj2,'PhasorSS') && isa(obj,'PhasorSS')
                if obj.isStatic 
                    % first check compatibility of the dimensions
                    if size(obj.D,2) ~= size(obj2.D,1)
                        error('The number of input of the first PhasorSS (size %d x %d) must be equal to the number of output of the second PhasorSS (size %d x %d)',size(obj.D,1),size(obj.D,2),size(obj2.D,1),size(obj2.D,2))
                    end

                    %first one is static, so we multiply the output of obj2 by the D matrix of obj
                    obj2.D = obj.D*obj2.D;
                    obj2.C = obj.D*obj2.C;
                    obj2.OutputName = obj.OutputName;
                    obj2.OutputUnit = obj.OutputUnit;
                    obj2.OutputGroup = obj.OutputGroup;
                    obj2 = obj2.checkIfReal();
                    objOut = obj2;
                    warning('p and T are propagated from the non-static PhasorSS, please check the values')
                    return
                end
                if obj2.isStatic
                    % first check compatibility of the dimensions
                    if size(obj.B,2) ~= size(obj2.D,1)
                        error('The number of input of the first PhasorSS must be equal to the number of output of the second PhasorSS')
                    end

                    %second one is static, so we multiply the input of obj by the D matrix of obj2
                    obj.B = obj.B*obj2.D;
                    obj.D = obj.D*obj2.D;
                    obj.InputName = obj2.InputName;
                    obj.InputUnit = obj2.InputUnit;
                    obj.InputGroup = obj2.InputGroup;
                    obj = obj.checkIfReal();
                    objOut = obj;
                    warning('p and T are propagated from the non-static PhasorSS, please check the values')
                    return
                end
                %else we are in full phasorSS case
                % check dimension compatibility
                if size(obj.B,2) ~= size(obj2.C,1)
                    error('The number of inputs of the first PhasorSS (%d) must be equal to the number of outputs of the second PhasorSS (%d)',size(obj.B,2),size(obj2.C,1))
                end

                newA = [obj.A obj.B*obj2.C;...
                 zeros(size(obj2.A,1),size(obj.A,2)) obj2.A];
                newB = [obj.B*obj2.D;obj2.B];
                newC = [obj2.C obj.D*obj2.C];
                newD = obj.D*obj2.D;

                %so obj2 output are the input of obj
                %obj output are the output of new system
                %obj2 input are the input of new system

                objOut = PhasorSS(newA,newB,newC,newD,obj.T,...
                    'StateName',[obj.StateName obj2.StateName],...
                    'StateUnit',[obj.StateUnit obj2.StateUnit],...
                    'InputName',obj2.InputName,...
                    'InputUnit',obj2.InputUnit,...
                    'InputGroup',obj2.InputGroup,...
                    'OutputName',obj.OutputName,...
                    'OutputUnit',obj.OutputUnit,...
                    'OutputGroup',obj.OutputGroup); 
                % 'Name',obj.Name,... 'Notes',obj.Notes,... 'UserData',obj.UserData,... 'p',obj.p were not empty in obj or obj2, issue a warning
                if ~isempty(obj.Name) || ~isempty(obj.Notes) || ~isempty(obj.UserData) || ~isempty(obj.p) || ~isempty(obj2.Name) || ~isempty(obj2.Notes) || ~isempty(obj2.UserData) || ~isempty(obj2.p)
                    warning('Name, Notes, UserData and p are not propagated from the PhasorSS objects, please check the values')
                end

                %if one of obj.p or obj2.p is not empty and the other is empty, issue a warning and propagate the non empty one
                if xor(isempty(obj.p),isempty(obj2.p))
                    if isempty(obj.p)
                        objOut.p = @(t,x,u) obj2.p(t,x((size(obj.A,1)+1):end),u);
                        objOut.isLPV = obj2.isLPV;
                        checkp(objOut);
                        warning('p is propagated from right PhasorSS, please check the values')
                    else
                        warning('Cannot propagate p from left PhasorSS, please check the values and update p manually')
                    end
                end

            end

        end

        function out = feedback(obj,obj2,varg)
            %FEEDBACK Feedback connection of two PhasorSS objects
            %   out = FEEDBACK(obj,obj2) computes the negative output feedback connection of two PhasorSS objects.
            %
            %
            % Negative Feedback Connection:
            %
            %      +      +------+
            %    u--o --->| obj  |-------> y
            %       |-    +------+    |
            %       |                 |
            %       |     +------+    |
            %       +-----| obj2 |<---+
            %             +------+
            %
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       obj2 - Instance of the PhasorSS class, or PhasorArray class, or matrix, if empty, obj2 is set to 1 (unary negative feedback)
            %   Optionnal name-value pair:
            %       'feedbackInput' - Feedback input, either 'state' or 'output' (default is 'output'). If 'state', the feedback is taken on the state, if 'output', the feedback is taken on the output.
            %       'uIndex' - Index of the input of obj to be used as feedback input (default is [])
            %       'yIndex' - Index of the output of obj2 to be used as feedback input (default is [])
            %       if "state" is specified, "yIndex" is ignored, a warning is issued if "yIndex" is not empty
            %
            %   Outputs:
            %       out - PhasorSS object
            %

            %handle case where obj2 is simply a matrix, convert it to static PhasorSS
            arguments
                obj
                obj2 = 1;
                varg.feedbackInput {mustBeMember(varg.feedbackInput,{'state','output'})} = 'output'
                varg.uIndex = []
                varg.yIndex = []
            end

            if isa(obj2, 'double') || isa(obj2,'PhasorArray')
                if isscalar(obj2)
                    obj2 = eye(obj.nu)*obj2;
                end
                obj2 = PhasorSS([], [], [], obj2,obj.T);
            end

            switch varg.feedbackInput
                case 'state'
                    if ~isempty(varg.yIndex)
                        warning('yIndex is ignored when feedbackInput is set to "state"')
                    end

                    %create a copy of phasorSS from obj where the output is the state instead of OG output
                    objState = PhasorSS(obj.A,obj.B,...
                        [],[],...
                        obj.T,...
                        'StateName',obj.StateName,...
                        'StateUnit',obj.StateUnit,...
                        'InputName',obj.InputName,...
                        'InputUnit',obj.InputUnit,...
                        'InputGroup',obj.InputGroup,...
                        'OutputName',obj.StateName,...
                        'OutputUnit',obj.StateUnit,...
                        'OutputGroup',obj.OutputGroup);
                        objState.p = obj.p;
                        objState.isLPV = obj.isLPV;


                    %perform the output feedback between objState and obj2
                    out = feedback(objState,obj2,"feedbackInput","output");

                    %retrieve the original output
                    newC = [obj.C-obj.D*obj2.D -obj.D*obj2.C];
                    newD = [obj.D];

                    out.C = newC;
                    out.D = newD;
                    out.OutputName = obj.OutputName;
                    out.OutputUnit = obj.OutputUnit;
                    out.OutputGroup = obj.OutputGroup;
                    return
                case 'output'
                    %ok je crois
            end


            if obj.T ~= obj2.T
                error('The periods of the two PhasorSS objects must be equal')
            end


            A1 = obj.A;
            B1 = obj.B;
            C1 = obj.C;
            D1 = obj.D;
            A2 = obj2.A;
            B2 = obj2.B;
            C2 = obj2.C;
            D2 = obj2.D;

            %check dimension compatibility
            if size(D1,2) ~= size(D2,1)
                error('The number of input of the first PhasorSS (size %d x %d) must be equal to the number of output of the second PhasorSS (size %d x %d)',size(D1,1),size(D1,2),size(D2,1),size(D2,2))
            end
            if size(D1,1) ~= size(D2,2)
                error('The number of output of the first PhasorSS (size %d x %d) must be equal to the number of input of the second PhasorSS (size %d x %d)',size(D1,1),size(D1,2),size(D2,1),size(D2,2))
            end

            %issue a warning about inversion if D1D2 is nonzero
            

            %compute the inverse
            invD1D2 = reduce(inv(eye(size(D1,1)) + D1*D2));


            newA = [A1-B1*D2*invD1D2*C1 -B1*C2+B1*D2*invD1D2*D1*C2;...
                B2*invD1D2*C1 A2-B2*invD1D2*D1*C2];
            newB = [B1-B1*D2*invD1D2*D1;B2*invD1D2*D1];
            newC = [invD1D2*C1 -invD1D2*D1*C2];
            newD = invD1D2*D1;

            newA = reduce(newA);
            newB = reduce(newB);
            newC = reduce(newC);
            newD = reduce(newD);

            out = PhasorSS(newA,newB,newC,newD,obj.T,...
                'StateName',[obj.StateName obj2.StateName],...
                'StateUnit',[obj.StateUnit obj2.StateUnit],...
                'InputName',obj.InputName,...
                'InputUnit',obj.InputUnit,...
                'InputGroup',obj.InputGroup,...
                'OutputName',obj.OutputName,...
                'OutputUnit',obj.OutputUnit,...
                'OutputGroup',obj.OutputGroup);

        end

        function lft(sys1,sys2,nu,ny)
            %LFT Lower LFT or Upper LFT of two PhasorSS objects, according to redheffer star product
            %   sys = LFT(sys1,sys2) computes the lower LFT of two PhasorSS objects.
            %   sys = LFT(sys1,sys2,nu,ny) computes the lower LFT of two PhasorSS objects with nu inputs and ny outputs.
            %
            %
            %This feedback loop connects the first nu outputs of sys2 to the last nu inputs of sys1 (signals u),
            %  and the last ny outputs of sys1 to the first ny inputs of sys2 (signals y).
            % The resulting system sys maps the input vector [w1 ; w2] to the output vector [z1 ; z2].
            %
            % In sys1, input is partionned as follows [w1 ; u1] (last nu inputs are u1, the control inputs)
            %        , output is partionned as follows [z1 ; y1] (last ny outputs are y1, the mesured outputs)
            % In sys2, input is partionned as follows [u2 ; w2] (first nu inputs are u2, the control inputs)
            %        , output is partionned as follows [y2 ; z2] (first ny outputs are y2, the mesured outputs)
            %
            %   Effectively, nu = dim(u1) = dim(y2) and ny = dim(y1) = dim(u2).
            %
            %  The resulting system sys maps the input vector [w1 ; w2] to the output vector [z1 ; z2].
            %
            %   The feedback loop is as follows:
            %                  +------+               
            %      w1 -------->| sys1 |------------> z1    
            %                  |      |           
            %           u1 --->|      |-------> y1
            %           |      +------+        | 
            %           |                      |   
            %           |      +------+        |  
            %           y2 <---| sys2 |<------ u2
            %                  |      |            
            %      z2 <--------|      |<------------ w2
            %                  +------+           
            %
            % sys = lft(sys1,sys2) produces:
            % The lower LFT of sys1 and sys2 if sys2 has fewer inputs and outputs than sys1.
            %   This amounts to deleting w2 and z2 in the above diagram.
            %
            % The upper LFT of sys1 and sys2 if sys1 has fewer inputs and outputs than sys2.
            %   This amounts to deleting w1 and z1 in the above diagram.
            
                nu1_base = sys1.nu;
                ny1_base = sys1.ny;
                nu2_base = sys2.nu;
                ny2_base = sys2.ny;
                if nargin < 3                    
                    if nu2_base<=ny1_base && ny2_base<=nu1_base
                        %no z and w for sys 2, lower lft, sys2 determine the number of inputs and outputs
                        nu = ny2_base;
                        ny = nu2_base;
                    elseif nu1_base<=ny2_base && ny1_base<=nu2_base
                        %no z and w for sys 1, upper lft, sys1 determine the number of inputs and outputs
                        nu = nu1_base;
                        ny = ny1_base;
                    else
                        error('The number of inputs and outputs of sys1 and sys2 must be such that one of them has fewer inputs and outputs than the other, or specify nu and ny')
                    end
                elseif nargin == 3
                    error('Specify both nu and ny')
                end                

                ny1 = ny;
                nz1 = ny1_base - ny1;
                nu1 = nu;
                nw1 = nu1_base - nu1;

                ny2 = nu1;
                nz2 = ny2_base - ny2;
                nu2 = ny1;
                nw2 = nu2_base - nu2;

                %check all are non negative
                if any([ny1 nz1 nu1 nw1 ny2 nz2 nu2 nw2]<0)
                    error('Check the values of nu and ny, negative value encountered')
                end

                %partitions the system as :  Part = [A B1 B2; C1 D11 D12; C2 D21 D22] 
                %                            [dx2;y2;z2] = Part2*[x2;u2;w2]
                %                            [dx1;z1;y1] = Part1*[x1;w1;u1]

                part1.A  = sys1.A;
                part1.B1 = sys1.B{:,1:nw1};
                part1.B2 = sys1.B{:,nw1+1:end};
                part1.C1 = sys1.C{1:nz1,:};
                part1.C2 = sys1.C{nz11+1:end,:};
                part1.D11 = sys1.D{1:nz1,1:nw1};
                part1.D12 = sys1.D{1:nz1,nw1+1:end};
                part1.D21 = sys1.D{ny1+1:end,1:nw1};
                part1.D22 = sys1.D{ny1+1:end,nw1+1:end};

                part2.A  = sys2.A;
                part2.B1 = sys2.B{:,1:nu2};
                part2.B2 = sys2.B{:,nu2+1:end};
                part2.C1 = sys2.C{1:ny2,:};
                part2.C2 = sys2.C{ny2+1:end,:};
                part2.D11 = sys2.D{1:ny2,1:nu2};
                part2.D12 = sys2.D{1:ny2,nu2+1:end};
                part2.D21 = sys2.D{ny2+1:end,1:nu2};
                part2.D22 = sys2.D{ny2+1:end,nu2+1:end};

                delta1 = -part1.D22*part2.D11+eye(ny1);
                delta2 = -part2.D11*part1.D22+eye(ny2);

                %compute the inverse
                if h(neglect(delta1))>0
                    warning('Δ₁ is periodic, inversion is involved, the result may be inaccurate')
                end
                if h(neglect(delta2))>0
                    warning('Δ₂ is periodic, inversion is involved, the result may be inaccurate')
                end
                invDelta1 = reduce(inv(neglect(delta1)));
                invDelta2 = reduce(inv(neglect(delta2)));

                %compute the new system
                newA11 = part1.A+part1.B2*invDelta2*part2.D11*part1.C2;
                newA12 = part1.B2*invDelta2*part2.C1;
                newA21 = part2.B1*invDelta1*part1.C2;
                newA22 = part2.A+part2.B1*invDelta1*part1.D22*part2.C1;
                newA = [newA11 newA12;newA21 newA22];

                newB11 = part1.B1+part1.B2*invDelta2*part2.D11*part1.D21;
                newB12 = part1.B2*invDelta2*part2.D12;
                newB21 = part2.b1*invDelta1*part1.D21;
                newB22 = part2.B2+part2.B1*invDelta1*part1.D22*part2.D12;
                newB = [newB11 newB12;newB21 newB22];

                newC11 = part1.C1 + part1.D12*invDelta2*part2.D11*part1.C2;
                newC12 = part1.D12*invDelta2*part2.C1;
                newC21 = part2.D21*invDelta1*part1.C2;
                newC22 = part2.C2 + part2.D21*invDelta1*part1.D22*part2.C1;
                newC = [newC11 newC12;newC21 newC22];

                newD11 = part1.D11 + part1.D12*invDelta2*part2.D11*part1.D21;
                newD12 = part1.D12*invDelta2*part2.D12;
                newD21 = part2.D21*invDelta1*part1.D21;
                newD22 = part2.D22 + part2.D21*invDelta1*part1.D22*part2.D12;
                newD = [newD11 newD12;newD21 newD22];

                w1InputName = sys1.InputName(1:nw1);
                w2InputName = sys2.InputName(nu2+1:end);
                newInputName = [w1InputName w2InputName];
                newInputUnit = [sys1.InputUnit(1:nw1) sys2.InputUnit(nu2+1:end)];


                z1OutputName = sys1.OutputName(1:nz1);
                z2OutputName = sys2.OutputName(ny2+1:end);
                newOutputName = [z1OutputName z2OutputName];
                newOutputUnit = [sys1.OutputUnit(1:nz1) sys2.OutputUnit(ny2+1:end)];

                out = PhasorSS(newA,newB,newC,newD,sys1.T,...
                    'StateName',[sys1.StateName sys2.StateName],...
                    'StateUnit',[sys1.StateUnit sys2.StateUnit],...
                    'InputName',newInputName,...
                    'InputUnit',newInputUnit,...
                    'OutputName',newOutputName,...
                    'OutputUnit',newOutputUnit);



        end
        
        function [A,B,C,D] = evalAngle(phasorSSobject,angle)
            A = phasorSSobject.A;
            B = phasorSSobject.B;
            C = phasorSSobject.C;
            D = phasorSSobject.D;
            if phasorSSobject.isReal
                A = A.SinCosForm();
                B = B.SinCosForm();
                C = C.SinCosForm();
                D = D.SinCosForm();
                h = (size(A,3)-1)/2;
                eit=[sin((h:-1:1)'*angle); cos((0:h)'*angle) ];
                A=tensorprod(A,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
                B=tensorprod(B,double(eit),3,1);
                C=tensorprod(C,double(eit),3,1);
                D=tensorprod(D,double(eit),3,1);
            else
                h = (size(Phas,3)-1)/2;
                eit=exp(1i*(-h:h)'*angle);
                A=tensorprod(A,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
                B=tensorprod(B,double(eit),3,1);
                C=tensorprod(C,double(eit),3,1);
                D=tensorprod(D,double(eit),3,1);
            end
        end

        function stem(phasorSSobject)
            T = tiledlayout(2,2);
            stem(phasorSSobject.A,'parent',T)
            title('A')
            stem(phasorSSobject.B,'parent',T)
            title('B')
            stem(phasorSSobject.C,'parent',T)
            title('C')
            stem(phasorSSobject.D,'parent',T)
            title('D')
        end

        function out = trunc(obj,h)
            %TRUNC Truncate the PhasorSS object to h harmonics
            %   TRUNC(obj,h) truncates the PhasorSS object to h harmonics, ie each matrix A,B,C,D
            %   is truncated to h harmonics.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       h - Number of harmonics to consider (truncation) (integer)
            %
            %   Outputs:
            %       out - Truncated PhasorSS object
            %

            out = obj;
            out.A = obj.A.trunc(h);
            out.B = obj.B.trunc(h);
            out.C = obj.C.trunc(h);
            out.D = obj.D.trunc(h);
        end

        function out = neglect(obj,threshold)
            %NEGLECT Neglect the PhasorSS object below a threshold
            %   NEGLECT(obj,threshold) neglects the PhasorSS object below a threshold, ie in each matrix A,B,C,D
            %   phasors below the threshold are set to zero.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %       threshold - Threshold to neglect (double)
            %
            %   Outputs:
            %       out - Neglected PhasorSS object
            %

            out = obj;
            out.A = obj.A.neglect(threshold);
            out.B = obj.B.neglect(threshold);
            out.C = obj.C.neglect(threshold);
            out.D = obj.D.neglect(threshold);
        end

        function plot(obj)
            %PLOT Plot the PhasorSS object
            %   PLOT(obj) plots the PhasorSS object calling the plot method of PhasorArray class for each matrix A,B,C,D.
            %

            %Concate the matrices A,B,C,D in a PhasorArray
            superMatrix = [obj.A , obj.B ; obj.C , obj.D];

            ff =gcf;
            % ff.Visible = 'off';
            if obj.isLPV
                T = 2*pi;
            else
                T= obj.T;
            end
            superMatrix.plot(T);
            n1 = obj.nx + obj.nu;
            n2 = obj.ny + obj.nx;

            %find all children that are axes
            allAxes = findall(ff.Children,'type','Axes');

            %delete all YLabel and title of allAxes
            arrayfun(@(ax) delete([ax.YLabel, ax.Title]), allAxes);

            for ii = 1:(obj.nx)
                nexttile(ii)
                title(obj.StateName{ii})
                nexttile((ii-1)*n1+1)
                ylabel(obj.StateName{ii},"FontWeight","bold","Rotation",0)
            end
            for ii = 1:(obj.nu)
                nexttile(obj.nx+ii)
                title(obj.InputName{ii})
            end
            for ii = 1:obj.ny
                nexttile((obj.nx+ii-1)*n1+1)
                ylabel(obj.OutputName{ii},"FontWeight","bold","Rotation",0)
            end

            % allAxes = findall(ff.Children,'type','Axes');
            linkaxes(allAxes,'x');
            %linkaxes y by row
            for ii = 1:n2
                linkaxes(allAxes((ii-1)*n1+1:ii*n1),'y');
            end

            % ff.Visible = 'on';

        end 


    end

    methods (Access=protected)
        function header = getHeader(obj)
              header = sprintf('PhasorSS object with properties: \n');
        end
        
        function propgrp = getPropertyGroups(obj)
            if obj.isLPV
                propList = struct(  'A',obj.A,...
                                    'B',obj.B,...
                                    'C',obj.C,...
                                    'D',obj.D,...
                                    'StateName',{obj.StateName},...
                                    'InputName',{obj.InputName},...
                                    'OutputName',{obj.OutputName},...
                                    'p',obj.p);
            else
                propList = struct('A',obj.A,...
                'B',obj.B,...
                'C',obj.C,...
                'D',obj.D,...
                'StateName',{obj.StateName},...
                'InputName',{obj.InputName},...
                'OutputName',{obj.OutputName},...
                'T',obj.T);
            end
            propgrp = matlab.mixin.util.PropertyGroup(propList);
        end

        function footer = getFooter(obj)
            varName = inputname(1);
            % Create the footer with a hyperlink to display all properties
            if obj.isLPV
                footer = sprintf('Continuous time Linear Parameter-Varying State-Space model with %d states, %d inputs, %d outputs and parameter function p(t)', obj.nx, obj.nu, obj.ny);
            else
                footer = sprintf('Continuous time Linear Time-Varying State-Space model with %d states, %d inputs, %d outputs and period T=%g\n See <a href="matlab:details(%s)">all properties</a>', obj.nx, obj.nu, obj.ny, obj.T, varName);
            end
        end
        
        

        function result = parenReference(obj, indexOp)
            % Handle A(I, J) indexing
            if numel(indexOp) == 1 && strcmp(indexOp(1).Type, 'Paren')
                subs = indexOp(1).Indices;
                if length(subs) == 2
                    I = subs{1}; % Output indices
                    J = subs{2}; % Input indices

                    % Extract submatrices
                    B_sub = obj.B{:,J};
                    C_sub = obj.C{I,:};
                    D_sub = obj.D{I,J};

                    % Handle empty cases for properties
                    InputName = obj.InputName;
                    if isempty(InputName)
                        InputName = repmat({''}, 1, size(B_sub, 2));
                    else
                        InputName = InputName(J);
                    end

                    InputUnit = obj.InputUnit;
                    if isempty(InputUnit)
                        InputUnit = repmat({''}, 1, size(B_sub, 2));
                    else
                        InputUnit = InputUnit(J);
                    end

                    OutputName = obj.OutputName;
                    if isempty(OutputName)
                        OutputName = repmat({''}, 1, size(C_sub, 1));
                    else
                        OutputName = OutputName(I);
                    end

                    OutputUnit = obj.OutputUnit;
                    if isempty(OutputUnit)
                        OutputUnit = repmat({''}, 1, size(C_sub, 1));
                    else
                        OutputUnit = OutputUnit(I);
                    end

                    % Create a new PhasorSS instance with the submatrices
                    result = PhasorSS(obj.A, B_sub, C_sub, D_sub, obj.T, ...
                        'isReal', obj.isReal, ...
                        'StateName', obj.StateName, ...
                        'StateUnit', obj.StateUnit, ...
                        'InputName', InputName, ...
                        'InputUnit', InputUnit, ...
                        'InputGroup', obj.InputGroup, ...
                        'OutputName', OutputName, ...
                        'OutputUnit', OutputUnit, ...
                        'OutputGroup', obj.OutputGroup, ...
                        'Name', obj.Name, ...
                        'Notes', obj.Notes, ...
                        'UserData', obj.UserData, ...
                        'p', obj.p);
                else
                    error('Invalid indexing. Use A(I, J) to extract a sub-system.');
                end
            else
                % Default behavior for other types of indexing
                result = builtin('subsref', obj, indexOp);
            end
        end

        function n = parenListLength(obj, indexOp, indexingContext)
            % Return the number of elements in the result of the indexing operation
            if numel(indexOp) == 1 && strcmp(indexOp(1).Type, '()')
                subs = indexOp(1).Indices;
                if length(subs) == 2
                    n = 1; % Always return 1 for A(I, J) indexing
                else
                    error('Invalid indexing. Use A(I, J) to extract a sub-system.');
                end
            else
                % Default behavior for other types of indexing
                n = builtin('numel', obj);
            end
        end

        function parenAssign(obj,indexOp,value)
            error('Parentheses assignment is not supported for PhasorSS objects. Use dot notation to assign properties.');
        end

        function parenDelete(obj,indexOp)
            error('Parentheses deletion is not supported for PhasorSS objects. Use dot notation to delete properties.');
        end
    end
    methods (Access=public)
        function varargout = size(obj,index)
            %SIZE Size of the PhasorSS object
            %   S = SIZE(obj) returns the size of the PhasorSS object.
            %
            %   Inputs:
            %       obj - Instance of the PhasorSS class
            %
            %   Outputs:
            %       S - Size of the PhasorSS object
            %   Example:
            %       S = size(obj);
            %       Return the size of the PhasorSS object
            %


            if nargin == 1
                if nargout ==0
                    fprintf('Periodic State-Space model with %d outputs, %d inputs and %d states\n',size(obj.C,1),size(obj.B,2),size(obj.A,1));
                else
                varargout{1} = [size(obj.C,1),size(obj.B,2)];
                end
            else
                if index == 1
                    varargout{1} = size(obj.C,1);
                elseif index == 2
                    varargout{1} = size(obj.B,2);
                else
                    error('Index out of range');
                end
            end
        end

        function obj = cat(dim,varargin)
            error('Concatenation is not supported for PhasorSS objects.');
        end
    end
    
    methods (Static, Access=public)
        function obj = empty()
            %EMPTY Create an empty PhasorSS object
            %   obj = EMPTY() creates an empty PhasorSS object.
            %
            %   Outputs:
            %       obj - Empty instance of the PhasorSS class
            %
            %   Example:
            %       obj = empty();
            %       Create an empty PhasorSS object
            %
            %   See also: ZEROS, ONES
            obj = PhasorSS([],[],[],[],[],'isReal',true);
        end

        function obj = fromSS(sys,T)
            %FROMSS Create a PhasorSS object from a State-Space object
            %   obj = FROMSS(sys,T) creates a PhasorSS object from a State-Space object.
            %
            %   Inputs:
            %       sys - State-Space object
            %       T - Period of the PhasorSS object
            %
            %   Outputs:
            %       obj - Instance of the PhasorSS class, with the same properties as the State-Space object
            arguments
                sys
                T = []
            end

            %if sys is a PhasorSS object, return it
            if isa(sys,'PhasorSS')
                obj = sys;
                return
            end
            %if sys is a transfer function convert it to a state space
            if isa(sys,'tf') 
                sys = ss(sys);
            end
            
            A = sys.A;
            B = sys.B;
            C = sys.C;
            D = sys.D;
            
            if isempty(T)
            T = 1;
            warning('Using FROMSS : Period T is empty, set a value to use the system as LTV')
            end
            
            obj = PhasorSS(A, B, C, D, T, 'isReal', isreal(sys), ...
            'StateName', sys.StateName, 'StateUnit', sys.StateUnit, ...
            'InputName', sys.InputName, 'InputUnit', sys.InputUnit, 'InputGroup', sys.InputGroup, ...
            'OutputName', sys.OutputName, 'OutputUnit', sys.OutputUnit, 'OutputGroup', sys.OutputGroup, ...
            'Name', sys.Name, 'Notes', sys.Notes, 'UserData', sys.UserData);
        end
    end
    
end


function Mt = phasor2time(Phas,angle)
    %PHASOR2TIME evaluate the PhasorArray in the time domain considering a phase angle and a complex valued phasorArray
    %   Mt = PHASOR2TIME(Phas,angle) evaluates the PhasorArray in the time domain considering a phase angle.
    %
    %   Inputs:
    %       Phas - PhasorArray
    %       angle - Phase angle
    %
    %   Outputs:
    %       Mt - Time domain evaluation of the PhasorArray
    %
    %   Example:
    %       M = PhasorArray(cat(3,ones(2,2),zeros(2,2),ones(2,2))); (2 cos(2 θ) + 2 cos(θ) + 0)
    %       angle = pi/2;
    %       Mt = phasor2time(M,angle);
    %       Evaluate the PhasorArray in the time domain considering a phase angle
    %
    %   See also: sincos2time
    h = (size(Phas,3)-1)/2;
    eit=exp(1i*(-h:h)'*angle);
    try
        Mt=tensorprod(Phas,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
    catch
        Mt=0;
        for ii = 1:size(Phas,3)
            Mt = Mt + Phas(:,:,ii)*eit(ii);
        end
    end
end


function Mt = sincos2time(PhasSC,angle)
    %SINCOS2TIME evaluate the PhasorArray in the time domain considering a phase angle and a SinCos valued phasorArray
    %   Mt = SINCOS2TIME(PhasSC,angle) evaluates the PhasorArray in the time domain considering a phase angle.
    %
    %   Inputs:
    %       PhasSC - SinCos PhasorArray
    %       angle - Phase angle
    %
    %   Outputs:
    %       Mt - Time domain evaluation of the PhasorArray
    %
    %   Example:
    %       M = PhasorArray(cat(3,ones(2,2),zeros(2,2),ones(2,2))); (2 cos(2 θ) + 2 cos(θ) + 0)
    %       angle = pi/2;
    %       MSC = M.SinCosForm();
    %       Mt = sincos2time(MSC,angle);
    %       Evaluate the PhasorArray in the time domain considering a phase angle
    %
    %   See also: phasor2time

    h = (size(PhasSC,3)-1)/2;
    eit=[sin((h:-1:1)'*angle); cos((0:h)'*angle) ];
    Mt=tensorprod(PhasSC,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
end

function structout = mergeStruct(struct1,struct2)
    %merge struct helps to merge 2 outputGroups structure, for exemple when are addind outputs to an already existing phasorSS
    %   mergeStruct(struct1,struct2) merge the struct2 into the struct1
    %
    %   Inputs:
    %       struct1 - First structure
    %       struct2 - Second structure
    %
    %   Outputs:
    %       structout - Merged structure 

    %merge the 2 structures
    structout = struct1;
    fn = fieldnames(struct2);
    for i = 1:numel(fn)
        if isfield(struct1,fn{i})
            structout.(fn{i}) = [struct1.(fn{i}) struct2.(fn{i})];
        else
            structout.(fn{i}) = struct2.(fn{i});
        end
    end
end