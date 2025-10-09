classdef sparsePhasorArray < PhasorArray
    %SPARSEPHASORARRAY Class for handling sparse phasor arrays
    %   This class inherits from PhasorArray and adds functionality for
    %   handling sparse phasor arrays.

    properties
        phasorIndex % Index indicating the phasor order for each page of Phasor3D
    end

    methods
        function obj = sparsePhasorArray(Phasor3D, phasorIndex)
            % Constructor for sparsePhasorArray
            %   obj = sparsePhasorArray(Phasor3D, phasorIndex) creates a
            %   sparsePhasorArray from a 3D array Phasor3D and a vector
            %   phasorIndex indicating the phasor order for each page.
            obj@PhasorArray(Phasor3D);
            obj.phasorIndex = phasorIndex;
        end

        function PA = toPhasorArray(obj)
            % Convert sparsePhasorArray to PhasorArray
            %   PA = toPhasorArray(obj) converts the sparsePhasorArray to
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

        function varargout = subsref(obj, S)
            % Override subsref to handle sparsePhasorArray
            %   This method calls the corresponding PhasorArray method
            %   after converting the sparsePhasorArray to a PhasorArray.
            PA = obj.toPhasorArray();
            [varargout{1:nargout}] = subsref(PA, S);
        end

        function varargout = plus(obj, other)
            % Override plus to handle sparsePhasorArray
            %   This method calls the corresponding PhasorArray method
            %   after converting the sparsePhasorArray to a PhasorArray.
            PA = obj.toPhasorArray();
            [varargout{1:nargout}] = plus(PA, other);
        end

        % Add similar overrides for other methods as needed
    end
end