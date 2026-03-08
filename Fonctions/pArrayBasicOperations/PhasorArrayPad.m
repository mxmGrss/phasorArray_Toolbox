function r = PhasorArrayPad(obj, delta_h)
    %PHASORARRAYPAD Pad a PhasorArray object with zeros in the harmonic domain.
    %   r = PhasorArrayPad(obj, delta_h) adds delta_h zero harmonics to both
    %   ends of the harmonic spectrum of the PhasorArray obj.
    %
    %   Input Arguments:
    %   - obj (PhasorArray) : The PhasorArray object to pad.
    %   - delta_h (integer) : If scalar, adds delta_h zeros at each end of the 
    %                         harmonic dimension (3rd dimension).
    %
    %   Output Arguments:
    %   - r (PhasorArray)   : The padded PhasorArray object.
    
    val = obj.value;
    [n1, n2, ~] = size(val);
    
    if isscalar(delta_h)
        if delta_h > 0
            % Add zeros to both ends of the 3rd dimension
            z = zeros(n1, n2, delta_h);
            val = cat(3, z, val, z);
        end
    elseif length(delta_h) == 3
        % Specialized padding if delta_h is a vector [h1 h2 h3]
        % (Legacy support based on PhasorArray.pad comments)
        if delta_h(1) > 0, val = cat(1, zeros(delta_h(1), n2, size(val,3)), val, zeros(delta_h(1), n2, size(val,3))); end
        if delta_h(2) > 0, val = cat(2, zeros(size(val,1), delta_h(2), size(val,3)), val, zeros(size(val,1), delta_h(2), size(val,3))); end
        if delta_h(3) > 0
            z = zeros(size(val,1), size(val,2), delta_h(3));
            val = cat(3, z, val, z);
        end
    end
    
    r = PhasorArray(val);
end
