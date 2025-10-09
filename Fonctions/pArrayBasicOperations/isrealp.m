function [r,R] = isrealp(A,tol)
    %ISREALP Checks if the imaginary part of the time-realization of a PhasorArray is negligible.
    %   [r, R] = isrealp(A, tol) checks if the imaginary part of the PhasorArray A is negligibly small
    %   compared to a given tolerance `tol`. It returns:
    %   - r: A boolean value indicating if A is real (true) or not (false).
    %   - R: A logical array showing element-wise comparison results.
    %
    %   Inputs:
    %       - A: A PhasorArray, symbolic or optimization variable to check.
    %       - tol (optional): Tolerance for comparison, default is an empty value, meaning it uses
    %         machine epsilon for comparison.
    %
    %   Outputs:
    %       - r: A boolean indicating if A is real.
    %       - R: A logical array showing element-wise comparison of real and imaginary parts.
    %
    %   See also: ismembertol, real, imag, flip, eps
arguments
    A
    tol=[]
end

if isa(A,"PhasorArray")
    A = A.value;
end

A_r  = real(A + flip(A,3))*(1/2) + 1i * imag(A - flip(A,3))*(1/2);

%if A is not sym, sdpvar nor ndsdpvar
if ~(isa(A,"sym") || isa(A,"ndsdpvar") || isa(A,"sdpvar"))
    if isempty(tol)
        r1 = ismembertol(real(A),real(A_r));
        r2 = ismembertol(imag(A),imag(A_r));
    else
        r1 = ismembertol(real(A),real(A_r),tol);
        r2 = ismembertol(imag(A),imag(A_r),tol);
    end
else %in this case ismembertol doesn't accept sym and sdpvar input, use manual difference instead
    if isempty(tol)
        r1 = abs(real(A)-real(A_r))<eps;
        r2 = abs(imag(A)-imag(A_r))<eps;
    else
        r1 = abs(real(A)-real(A_r))<tol;
        r2 = abs(imag(A)-imag(A_r))<tol;
    end
end
R = zeros(size(A),'logical');
R((r1 + r2)==2)= true;
r = all(R,'all');
end