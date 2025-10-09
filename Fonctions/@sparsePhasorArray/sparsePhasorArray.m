classdef sparsePhasorArray 
    %SPARSEPHASORARRAY Class for handling sparse phasor arrays
    %   This class inherits from PhasorArray and adds functionality for
    %   handling sparse phasor arrays.
    %   Every method from PhasorArray is available in this class through the 
    %   subsref method. If a method is not found in this class, the object is
    %   converted to a PhasorArray and the method is called on the PhasorArray
    %   object.
    %
    %   to call a method prefer the following syntax:
    %   result = obj.methodName(varargin)
    %
    %   the synthax methodName(obj,varargin) is not available for the methods from PhasorArray
    %   that are not overloaded in this class.

    properties
        phasorIndex % Index indicating the phasor order for each page of Phasor3D
        Phasor3D % 3D array of phasors
    end

    methods
        function obj = sparsePhasorArray(Phasor3D, phasorIndex)
            % Constructor for sparsePhasorArray
            %   obj = sparsePhasorArray(Phasor3D, phasorIndex) creates a
            %   sparsePhasorArray from a 3D array Phasor3D and a vector
            %   phasorIndex indicating the phasor order for each page.
            %
            %   In the case of a single input argument A, the constructor
            %   looks for zeros pages in A, deletes them, and determine the 
            %   phasorIndex accordingly.
            switch nargin
                case 1
                    if isa(Phasor3D, 'PhasorArray')
                        Phasor3D = Phasor3D.value;
                    end
                    % Find zero pages in Phasor3D
                    zeroPages = squeeze(all(all(Phasor3D == 0, 1), 2));
                    phasorIndex = find(~zeroPages) - ceil(length(zeroPages)/2);
                    % Remove zero pages
                    Phasor3D = Phasor3D(:, :, ~zeroPages);
                case 2
                    % Do nothing
                otherwise
                    error('Invalid number of input arguments');
            end
            obj.Phasor3D = Phasor3D;
            if iscolumn(phasorIndex)
                phasorIndex = phasorIndex';
            end
            obj.phasorIndex = phasorIndex;
        end

        function PA = toPhasorArray(obj)
            % Convert sparsePhasorArray to PhasorArray
            %   PA = PhasorArray(obj) converts the sparsePhasorArray to
            %   a regular PhasorArray by correctly assigning the pages
            %   according to the phasorIndex and filling zeros elsewhere.
            maxIndex = max(abs(obj.phasorIndex));
            n = size(obj.Phasor3D, 1);
            m = size(obj.Phasor3D, 2);
            Phasor3D_full = zeros(n, m, 2*maxIndex + 1);
            for i = 1:length(obj.phasorIndex)
                idx = obj.phasorIndex(i) + maxIndex + 1;
                Phasor3D_full(:, :, idx) = obj.Phasor3D(:, :, i);
            end
            PA = PhasorArray(Phasor3D_full);
        end

        function result = callMethod(obj, methodName, varargin)
            % Call a method on the PhasorArray after converting
            %   result = callMethod(obj, methodName, varargin) converts the
            %   sparsePhasorArray to a PhasorArray and calls the specified
            %   method with the given arguments.
            PA = obj.toPhasorArray();
            result = PA.(methodName)(varargin{:});
        end

        function out = subsref(obj, S)
            % Handle dot notation (spA_obj.methodName())
            if strcmp(S(1).type, '.')
                methodName = S(1).subs;
                
                if ismethod(obj, methodName)
                    out = builtin('subsref', obj, S);
                elseif ismethod(PhasorArray, methodName)
                    objA = obj.toPhasorArray();
                    out = subsref(objA, S);
                else
                    error("Method '%s' not found in sparsePhasorArray or PhasorArray.", methodName);
                end
            else
                out = builtin('subsref', obj, S);
            end
        end

        function out = isa(obj, className)
            disp('heho')
            % Pretend sparsePhasorArray is a PhasorArray to enable function-call syntax
            if strcmp(className, 'PhasorArray')
                out = true;
            else
                out = builtin('isa', obj, className);
            end
        end

        % Operator overloading
        function result = plus(obj, other)
            % Override plus to handle sparsePhasorArray
            %   This method calls the corresponding PhasorArray method
            %   after converting the sparsePhasorArray to a PhasorArray.
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA + other;
            result = sparsePhasorArray(result);
        end

        function result = minus(obj, other)
            % Override minus to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA - other;
            result = sparsePhasorArray(result);
        end

        function result = mtimes(obj, other)
            % Override mtimes to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA * other;
            result = sparsePhasorArray(result);
        end

        function result = rdivide(obj, other)
            % Override rdivide to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA / other;
        end

        function result = ldivide(obj, other)
            % Override ldivide to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA \ other;
        end

        function result = power(obj, other)
            % Override power to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            result = PA .^ other;
            result = sparsePhasorArray(result);
        end

        function result = mpower(obj, other)
            % Override mpower to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            result = PA ^ other;
            result = sparsePhasorArray(result);
        end

        function result = times(obj, other)
            % Override times to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA .* other;
            result = sparsePhasorArray(result);
        end

        function result = mrdivide(obj, other)
            % Override mrdivide to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA ./ other;
        end

        function result = mldivide(obj, other)
            % Override mldivide to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            if isa(other, 'sparsePhasorArray')
                other = other.toPhasorArray();
            end
            result = PA .\ other;
        end

        function result = uminus(obj)
            % Override uminus to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            result = -PA;
            result = sparsePhasorArray(result);
        end

        function result = uplus(obj)
            % Override uplus to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            result = +PA;
            result = sparsePhasorArray(result);
        end

        function result = transpose(obj)
            % Override transpose to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            result = PA.';
            result = sparsePhasorArray(result);
        end

        function result = ctranspose(obj)
            % Override ctranspose to handle sparsePhasorArray
            PA = obj.toPhasorArray();
            result = PA';
            result = sparsePhasorArray(result);
        end

        function out = getPhasorIndex(obj)
            % Get the phasor index
            out = obj.phasorIndex;
        end

        function val = evalp(obj, theta)
            % Evaluate the sparse phasor array at the given angles
            %   val = evalp(obj, theta) evaluates the sparse phasor array
            %   at the given angles theta.
            if iscolumn(theta)
                theta = theta';
            end
            if isreal(obj) 
                base = exp(1i * (obj.phasorIndex((end+1)/2:end)')* theta );
                val = real(tensorprod(cat(3,obj.Phasor3D(:,:,(end+1)/2),2*obj.Phasor3D(:,:,(1+(end+1)/2:end))), base,3,1));
            else
                base = exp(1i * (obj.phasorIndex')* theta );
                val = tensorprod(obj.Phasor3D, base,3,1);
            end
        end

        function out = isreal(obj,tol)
            arguments
                obj
                tol (1,1) double = 1e-12
            end
            % Check if the periodic matrix associated with the sparse phasor array is real

            objP = toPhasorArray(obj);
            out = objP.isreal(tol);
        end


    end

    methods (Static)
        function out = builtin(funcName, varargin)
            disp("heho")
            % Intercepts calls like methodName(spA_obj) if they fail
            try
                out = builtin(funcName, varargin{:});
            catch
                % Check if the first argument is sparsePhasorArray and the function is a PhasorArray method
                if numel(varargin) > 0 && isa(varargin{1}, 'sparsePhasorArray') && ismethod(PhasorArray, funcName)
                    objA = varargin{1}.toPhasorArray();
                    out = feval(funcName, objA, varargin{2:end});
                else
                    rethrow(lasterror());
                end
            end
        end
    end
end