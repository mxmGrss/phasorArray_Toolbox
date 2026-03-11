% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

nx = size(A,1);   nu = size(B,2);

% Speed selection row
C_omega = PhasorArray([0 0 1]);        % 1×3, picks ω_m

% Augmented system: x_tilde = [i_α; i_β; ω; z]
A_aug = [A,                      PhasorArray(zeros(nx,1));
         C_omega,                PhasorArray(0)];

B_aug = [B;
         PhasorArray(zeros(1,nu))];

nx_aug = size(A_aug, 1);   % = 4
