% +--------------------------------------------------------------------------+
% | 'TrussMe' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023 |
% | Current version authors:                                                 |
% |   Davide Stocco and Matteo Larcher.                                      |
% +--------------------------------------------------------------------------+

% Matlab generated code for system: Test_05
% This file has been automatically generated by TrussMe.
% DISCLAIMER: If you need to edit it, do it wisely!

classdef Test_05 < TrussMe.System
  %
  % Test_05 class
  %
  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function this = Test_05( varargin )
      % Class constructor.

      % User data
      if (nargin == 0)
        data = [];
      elseif (nargin == 1 && isstruct(varargin{1}))
        data = varargin{1};
      else
        error('wrong number of input arguments.');
      end

      % Call superclass constructor
      this = this@TrussMe.System(data);
    end % Test_05
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K = K( ~, ~, ~ )
      % Evaluate the stiffness matrix K.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_K = sparse(18, 18);
      out_K(1, 1) = 1000000000.;
      out_K(7, 1) = -1000000000.;
      out_K(2, 2) = 1500000.;
      out_K(6, 2) = 1500000.;
      out_K(8, 2) = -1500000.;
      out_K(12, 2) = 1500000.;
      out_K(3, 3) = 1500000.;
      out_K(5, 3) = -1500000.;
      out_K(9, 3) = -1500000.;
      out_K(11, 3) = -1500000.;
      out_K(4, 4) = 500000.;
      out_K(10, 4) = -500000.;
      out_K(3, 5) = -1500000.;
      out_K(5, 5) = 2000000.;
      out_K(9, 5) = 1500000.;
      out_K(11, 5) = 1000000.;
      out_K(2, 6) = 1500000.;
      out_K(6, 6) = 2000000.;
      out_K(8, 6) = -1500000.;
      out_K(12, 6) = 1000000.;
      out_K(1, 7) = -1000000000.;
      out_K(7, 7) = 2000000000.;
      out_K(13, 7) = -1000000000.;
      out_K(2, 8) = -1500000.;
      out_K(6, 8) = -1500000.;
      out_K(8, 8) = 3000000.;
      out_K(14, 8) = -1500000.;
      out_K(18, 8) = 1500000.;
      out_K(3, 9) = -1500000.;
      out_K(5, 9) = 1500000.;
      out_K(9, 9) = 3000000.;
      out_K(15, 9) = -1500000.;
      out_K(17, 9) = -1500000.;
      out_K(4, 10) = -500000.;
      out_K(10, 10) = 1000000.;
      out_K(16, 10) = -500000.;
      out_K(3, 11) = -1500000.;
      out_K(5, 11) = 1000000.;
      out_K(11, 11) = 4000000.;
      out_K(15, 11) = 1500000.;
      out_K(17, 11) = 1000000.;
      out_K(2, 12) = 1500000.;
      out_K(6, 12) = 1000000.;
      out_K(12, 12) = 4000000.;
      out_K(14, 12) = -1500000.;
      out_K(18, 12) = 1000000.;
      out_K(7, 13) = -1000000000.;
      out_K(13, 13) = 1000000000.;
      out_K(8, 14) = -1500000.;
      out_K(12, 14) = -1500000.;
      out_K(14, 14) = 1500000.;
      out_K(18, 14) = -1500000.;
      out_K(9, 15) = -1500000.;
      out_K(11, 15) = 1500000.;
      out_K(15, 15) = 1500000.;
      out_K(17, 15) = 1500000.;
      out_K(10, 16) = -500000.;
      out_K(16, 16) = 500000.;
      out_K(9, 17) = -1500000.;
      out_K(11, 17) = 1000000.;
      out_K(15, 17) = 1500000.;
      out_K(17, 17) = 2000000.;
      out_K(8, 18) = 1500000.;
      out_K(12, 18) = 1000000.;
      out_K(14, 18) = -1500000.;
      out_K(18, 18) = 2000000.;
    end % K
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_ff = K_ff( ~, ~, ~ )
      % Evaluate the stiffness matrix K_ff.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_K_ff = sparse(10, 10);
      out_K_ff(1, 1) = 2000000.;
      out_K_ff(4, 1) = 1000000.;
      out_K_ff(2, 2) = 1000000.;
      out_K_ff(8, 2) = -500000.;
      out_K_ff(3, 3) = 4000000.;
      out_K_ff(7, 3) = 1500000.;
      out_K_ff(9, 3) = 1000000.;
      out_K_ff(1, 4) = 1000000.;
      out_K_ff(4, 4) = 4000000.;
      out_K_ff(6, 4) = -1500000.;
      out_K_ff(10, 4) = 1000000.;
      out_K_ff(5, 5) = 1000000000.;
      out_K_ff(4, 6) = -1500000.;
      out_K_ff(6, 6) = 1500000.;
      out_K_ff(10, 6) = -1500000.;
      out_K_ff(3, 7) = 1500000.;
      out_K_ff(7, 7) = 1500000.;
      out_K_ff(9, 7) = 1500000.;
      out_K_ff(2, 8) = -500000.;
      out_K_ff(8, 8) = 500000.;
      out_K_ff(3, 9) = 1000000.;
      out_K_ff(7, 9) = 1500000.;
      out_K_ff(9, 9) = 2000000.;
      out_K_ff(4, 10) = 1000000.;
      out_K_ff(6, 10) = -1500000.;
      out_K_ff(10, 10) = 2000000.;
    end % K_ff
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_fs = K_fs( ~, ~, ~ )
      % Evaluate the stiffness matrix K_fs.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_K_fs = sparse(10, 8);
      out_K_fs(1, 2) = 1500000.;
      out_K_fs(4, 2) = 1500000.;
      out_K_fs(3, 3) = -1500000.;
      out_K_fs(2, 4) = -500000.;
      out_K_fs(3, 5) = 1000000.;
      out_K_fs(5, 6) = -1000000000.;
      out_K_fs(1, 7) = -1500000.;
      out_K_fs(6, 7) = -1500000.;
      out_K_fs(10, 7) = 1500000.;
      out_K_fs(7, 8) = -1500000.;
      out_K_fs(9, 8) = -1500000.;
    end % K_fs
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_sf = K_sf( ~, ~, ~ )
      % Evaluate the stiffness matrix K_sf.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_K_sf = sparse(8, 10);
      out_K_sf(2, 1) = 1500000.;
      out_K_sf(7, 1) = -1500000.;
      out_K_sf(4, 2) = -500000.;
      out_K_sf(3, 3) = -1500000.;
      out_K_sf(5, 3) = 1000000.;
      out_K_sf(2, 4) = 1500000.;
      out_K_sf(6, 5) = -1000000000.;
      out_K_sf(7, 6) = -1500000.;
      out_K_sf(8, 7) = -1500000.;
      out_K_sf(8, 9) = -1500000.;
      out_K_sf(7, 10) = 1500000.;
    end % K_sf
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_ss = K_ss( ~, ~, ~ )
      % Evaluate the stiffness matrix K_ss.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_K_ss = sparse(8, 8);
      out_K_ss(1, 1) = 1000000000.;
      out_K_ss(6, 1) = -1000000000.;
      out_K_ss(2, 2) = 1500000.;
      out_K_ss(7, 2) = -1500000.;
      out_K_ss(3, 3) = 1500000.;
      out_K_ss(5, 3) = -1500000.;
      out_K_ss(8, 3) = -1500000.;
      out_K_ss(4, 4) = 500000.;
      out_K_ss(3, 5) = -1500000.;
      out_K_ss(5, 5) = 2000000.;
      out_K_ss(8, 5) = 1500000.;
      out_K_ss(1, 6) = -1000000000.;
      out_K_ss(6, 6) = 2000000000.;
      out_K_ss(2, 7) = -1500000.;
      out_K_ss(7, 7) = 3000000.;
      out_K_ss(3, 8) = -1500000.;
      out_K_ss(5, 8) = 1500000.;
      out_K_ss(8, 8) = 3000000.;
    end % K_ss
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_d = d( ~, ~, ~ )
      % Evaluate the deformation vector d.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_d = sparse(18, 1);
      out_d(6) = .00333333333333333;
      out_d(12) = -.00666666666666666;
      out_d(14) = -.0266666666666666;
      out_d(18) = -.0166666666666666;
    end % d
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_d_f = d_f( ~, ~, ~ )
      % Evaluate the deformation vector d_f.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_d_f = sparse(10, 1);
      out_d_f(1) = .00333333333333333;
      out_d_f(4) = -.00666666666666666;
      out_d_f(6) = -.0266666666666666;
      out_d_f(10) = -.0166666666666666;
    end % d_f
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_d_s = d_s( ~, ~, ~ )
      % Evaluate the deformation vector d_s.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_d_s = sparse(8, 1);
    end % d_s
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f = f( ~, ~, ~ )
      % Evaluate the force vector f.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_f = sparse(18, 1);
      out_f(2) = -4999.99999999999;
      out_f(8) = 9999.99999999999;
      out_f(14) = -5000;
    end % f
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f_f = f_f( ~, ~, ~ )
      % Evaluate the force vector f_f.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_f_f = sparse(10, 1);
      out_f_f(6) = -5000;
    end % f_f
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f_s = f_s( ~, ~, ~ )
      % Evaluate the force vector f_s.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_f_s = sparse(8, 1);
      out_f_s(2) = -4999.99999999999;
      out_f_s(7) = 9999.99999999999;
    end % f_s
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f_r = f_r( ~, ~, ~ )
      % Evaluate the force vector f_r.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_f_r = sparse(8, 1);
    end % f_r
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_perm = perm( ~ )
      % Evaluate the permutation vector.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_perm = zeros(18, 1);
      out_perm(1) = 6;
      out_perm(2) = 10;
      out_perm(3) = 11;
      out_perm(4) = 12;
      out_perm(5) = 13;
      out_perm(6) = 14;
      out_perm(7) = 15;
      out_perm(8) = 16;
      out_perm(9) = 17;
      out_perm(10) = 18;
      out_perm(11) = 1;
      out_perm(12) = 2;
      out_perm(13) = 3;
      out_perm(14) = 4;
      out_perm(15) = 5;
      out_perm(16) = 7;
      out_perm(17) = 8;
      out_perm(18) = 9;
    end % perm
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_unperm = unperm( ~ )
      % Evaluate the unpermutation vector.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_unperm = zeros(18, 1);
      out_unperm(1) = 11;
      out_unperm(2) = 12;
      out_unperm(3) = 13;
      out_unperm(4) = 14;
      out_unperm(5) = 15;
      out_unperm(6) = 1;
      out_unperm(7) = 16;
      out_unperm(8) = 17;
      out_unperm(9) = 18;
      out_unperm(10) = 2;
      out_unperm(11) = 3;
      out_unperm(12) = 4;
      out_unperm(13) = 5;
      out_unperm(14) = 6;
      out_unperm(15) = 7;
      out_unperm(16) = 8;
      out_unperm(17) = 9;
      out_unperm(18) = 10;
    end % unperm
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_v = v( ~, ~ )
      % Evaluate the veiling vector.

      % Extract properties
      % No properties

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_v = zeros(0, 1);
    end % v
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
end % Test_05

% That's All Folks!
