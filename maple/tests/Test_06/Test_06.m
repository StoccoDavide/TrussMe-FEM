% +--------------------------------------------------------------------------+
% | 'TrussMe' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023 |
% | Current version authors:                                                 |
% |   Davide Stocco and Matteo Larcher.                                      |
% +--------------------------------------------------------------------------+

% Matlab generated code for system: Test_06
% This file has been automatically generated by TrussMe.
% DISCLAIMER: If you need to edit it, do it wisely!

classdef Test_06 < TrussMe.System
  %
  % Test_06 class
  %
  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function this = Test_06( varargin )
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
    end % Test_06
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
      out_K = sparse(24, 24);
      out_K(1, 1) = 201.3888889;
      out_K(7, 1) = -201.3888889;
      out_K(2, 2) = 87.40837192;
      out_K(6, 2) = 6293.402778;
      out_K(8, 2) = -87.40837192;
      out_K(12, 2) = 6293.402778;
      out_K(3, 3) = 87.40837192;
      out_K(5, 3) = -6293.402778;
      out_K(9, 3) = -87.40837192;
      out_K(11, 3) = -6293.402778;
      out_K(4, 4) = 151041.6667;
      out_K(10, 4) = -151041.6667;
      out_K(3, 5) = -6293.402778;
      out_K(5, 5) = 604166.6668;
      out_K(9, 5) = 6293.402778;
      out_K(11, 5) = 302083.3334;
      out_K(2, 6) = 6293.402778;
      out_K(6, 6) = 604166.6668;
      out_K(8, 6) = -6293.402778;
      out_K(12, 6) = 302083.3334;
      out_K(1, 7) = -201.3888889;
      out_K(7, 7) = 402.7777778;
      out_K(13, 7) = -201.3888889;
      out_K(2, 8) = -87.40837192;
      out_K(6, 8) = -6293.402778;
      out_K(8, 8) = 174.8167438;
      out_K(14, 8) = -87.40837192;
      out_K(18, 8) = 6293.402778;
      out_K(3, 9) = -87.40837192;
      out_K(5, 9) = 6293.402778;
      out_K(9, 9) = 174.8167438;
      out_K(15, 9) = -87.40837192;
      out_K(17, 9) = -6293.402778;
      out_K(4, 10) = -151041.6667;
      out_K(10, 10) = 302083.3334;
      out_K(16, 10) = -151041.6667;
      out_K(3, 11) = -6293.402778;
      out_K(5, 11) = 302083.3334;
      out_K(11, 11) = 1208333.334;
      out_K(15, 11) = 6293.402778;
      out_K(17, 11) = 302083.3334;
      out_K(2, 12) = 6293.402778;
      out_K(6, 12) = 302083.3334;
      out_K(12, 12) = 1208333.334;
      out_K(14, 12) = -6293.402778;
      out_K(18, 12) = 302083.3334;
      out_K(7, 13) = -201.3888889;
      out_K(13, 13) = 302.0833333;
      out_K(19, 13) = -100.6944444;
      out_K(8, 14) = -87.40837192;
      out_K(12, 14) = -6293.402778;
      out_K(14, 14) = 98.33441841;
      out_K(18, 14) = -4720.052084;
      out_K(20, 14) = -10.92604649;
      out_K(24, 14) = 1573.350694;
      out_K(9, 15) = -87.40837192;
      out_K(11, 15) = 6293.402778;
      out_K(15, 15) = 98.33441841;
      out_K(17, 15) = 4720.052084;
      out_K(21, 15) = -10.92604649;
      out_K(23, 15) = -1573.350694;
      out_K(10, 16) = -151041.6667;
      out_K(16, 16) = 226562.5000;
      out_K(22, 16) = -75520.83333;
      out_K(9, 17) = -6293.402778;
      out_K(11, 17) = 302083.3334;
      out_K(15, 17) = 4720.052084;
      out_K(17, 17) = 906250.0001;
      out_K(21, 17) = 1573.350694;
      out_K(23, 17) = 151041.6667;
      out_K(8, 18) = 6293.402778;
      out_K(12, 18) = 302083.3334;
      out_K(14, 18) = -4720.052084;
      out_K(18, 18) = 906250.0001;
      out_K(20, 18) = -1573.350694;
      out_K(24, 18) = 151041.6667;
      out_K(13, 19) = -100.6944444;
      out_K(19, 19) = 100.6944444;
      out_K(14, 20) = -10.92604649;
      out_K(18, 20) = -1573.350694;
      out_K(20, 20) = 10.92604649;
      out_K(24, 20) = -1573.350694;
      out_K(15, 21) = -10.92604649;
      out_K(17, 21) = 1573.350694;
      out_K(21, 21) = 10.92604649;
      out_K(23, 21) = 1573.350694;
      out_K(16, 22) = -75520.83333;
      out_K(22, 22) = 75520.83333;
      out_K(15, 23) = -1573.350694;
      out_K(17, 23) = 151041.6667;
      out_K(21, 23) = 1573.350694;
      out_K(23, 23) = 302083.3333;
      out_K(14, 24) = 1573.350694;
      out_K(18, 24) = 151041.6667;
      out_K(20, 24) = -1573.350694;
      out_K(24, 24) = 302083.3333;
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
      out_K_ff = sparse(15, 15);
      out_K_ff(1, 1) = 151041.6667;
      out_K_ff(7, 1) = -151041.6667;
      out_K_ff(2, 2) = 604166.6668;
      out_K_ff(6, 2) = 6293.402778;
      out_K_ff(8, 2) = 302083.3334;
      out_K_ff(3, 3) = 604166.6668;
      out_K_ff(5, 3) = -6293.402778;
      out_K_ff(9, 3) = 302083.3334;
      out_K_ff(4, 4) = 402.7777778;
      out_K_ff(3, 5) = -6293.402778;
      out_K_ff(5, 5) = 174.8167438;
      out_K_ff(12, 5) = 6293.402778;
      out_K_ff(2, 6) = 6293.402778;
      out_K_ff(6, 6) = 174.8167438;
      out_K_ff(11, 6) = -6293.402778;
      out_K_ff(1, 7) = -151041.6667;
      out_K_ff(7, 7) = 302083.3334;
      out_K_ff(10, 7) = -151041.6667;
      out_K_ff(2, 8) = 302083.3334;
      out_K_ff(8, 8) = 1208333.334;
      out_K_ff(11, 8) = 302083.3334;
      out_K_ff(3, 9) = 302083.3334;
      out_K_ff(9, 9) = 1208333.334;
      out_K_ff(12, 9) = 302083.3334;
      out_K_ff(7, 10) = -151041.6667;
      out_K_ff(10, 10) = 226562.5000;
      out_K_ff(13, 10) = -75520.83333;
      out_K_ff(6, 11) = -6293.402778;
      out_K_ff(8, 11) = 302083.3334;
      out_K_ff(11, 11) = 906250.0001;
      out_K_ff(14, 11) = 151041.6667;
      out_K_ff(5, 12) = 6293.402778;
      out_K_ff(9, 12) = 302083.3334;
      out_K_ff(12, 12) = 906250.0001;
      out_K_ff(15, 12) = 151041.6667;
      out_K_ff(10, 13) = -75520.83333;
      out_K_ff(13, 13) = 75520.83333;
      out_K_ff(11, 14) = 151041.6667;
      out_K_ff(14, 14) = 302083.3333;
      out_K_ff(12, 15) = 151041.6667;
      out_K_ff(15, 15) = 302083.3333;
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
      out_K_fs = sparse(15, 9);
      out_K_fs(4, 1) = -201.3888889;
      out_K_fs(3, 2) = 6293.402778;
      out_K_fs(5, 2) = -87.40837192;
      out_K_fs(9, 2) = 6293.402778;
      out_K_fs(2, 3) = -6293.402778;
      out_K_fs(6, 3) = -87.40837192;
      out_K_fs(8, 3) = -6293.402778;
      out_K_fs(4, 4) = -201.3888889;
      out_K_fs(5, 5) = -87.40837192;
      out_K_fs(9, 5) = -6293.402778;
      out_K_fs(12, 5) = -4720.052084;
      out_K_fs(15, 5) = 1573.350694;
      out_K_fs(6, 6) = -87.40837192;
      out_K_fs(8, 6) = 6293.402778;
      out_K_fs(11, 6) = 4720.052084;
      out_K_fs(14, 6) = -1573.350694;
      out_K_fs(12, 8) = -1573.350694;
      out_K_fs(15, 8) = -1573.350694;
      out_K_fs(11, 9) = 1573.350694;
      out_K_fs(14, 9) = 1573.350694;
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
      out_K_sf = sparse(9, 15);
      out_K_sf(3, 2) = -6293.402778;
      out_K_sf(2, 3) = 6293.402778;
      out_K_sf(1, 4) = -201.3888889;
      out_K_sf(4, 4) = -201.3888889;
      out_K_sf(2, 5) = -87.40837192;
      out_K_sf(5, 5) = -87.40837192;
      out_K_sf(3, 6) = -87.40837192;
      out_K_sf(6, 6) = -87.40837192;
      out_K_sf(3, 8) = -6293.402778;
      out_K_sf(6, 8) = 6293.402778;
      out_K_sf(2, 9) = 6293.402778;
      out_K_sf(5, 9) = -6293.402778;
      out_K_sf(6, 11) = 4720.052084;
      out_K_sf(9, 11) = 1573.350694;
      out_K_sf(5, 12) = -4720.052084;
      out_K_sf(8, 12) = -1573.350694;
      out_K_sf(6, 14) = -1573.350694;
      out_K_sf(9, 14) = 1573.350694;
      out_K_sf(5, 15) = 1573.350694;
      out_K_sf(8, 15) = -1573.350694;
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
      out_K_ss = sparse(9, 9);
      out_K_ss(1, 1) = 201.3888889;
      out_K_ss(2, 2) = 87.40837192;
      out_K_ss(3, 3) = 87.40837192;
      out_K_ss(4, 4) = 302.0833333;
      out_K_ss(7, 4) = -100.6944444;
      out_K_ss(5, 5) = 98.33441841;
      out_K_ss(8, 5) = -10.92604649;
      out_K_ss(6, 6) = 98.33441841;
      out_K_ss(9, 6) = -10.92604649;
      out_K_ss(4, 7) = -100.6944444;
      out_K_ss(7, 7) = 100.6944444;
      out_K_ss(5, 8) = -10.92604649;
      out_K_ss(8, 8) = 10.92604649;
      out_K_ss(6, 9) = -10.92604649;
      out_K_ss(9, 9) = 10.92604649;
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
      out_d = sparse(24, 1);
      out_d(6) = -.0113876724187004;
      out_d(8) = -1.36016586259839;
      out_d(12) = -.00556144396147923;
      out_d(14) = -1.5;
      out_d(18) = .00238344827777457;
      out_d(24) = .00662077585937338;
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
      out_d_f = sparse(15, 1);
      out_d_f(3) = -.0113876724187004;
      out_d_f(5) = -1.36016586259839;
      out_d_f(9) = -.00556144396147923;
      out_d_f(12) = .00238344827777457;
      out_d_f(15) = .00662077585937338;
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
      out_d_s = sparse(9, 1);
      out_d_s(5) = -1.5;
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
      out_f = sparse(24, 1);
      out_f(2) = 12.2222674792198;
      out_f(8) = -20;
      out_f(14) = 5.55546513530001;
      out_f(20) = 2.22226743988671;
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
      out_f_f = sparse(15, 1);
      out_f_f(5) = -20;
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
      out_f_s = sparse(9, 1);
      out_f_s(2) = 12.2222674792198;
      out_f_s(5) = 5.55546513530001;
      out_f_s(8) = 2.22226743988671;
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
      out_f_r = sparse(9, 1);
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
      out_perm = zeros(24, 1);
      out_perm(1) = 4;
      out_perm(2) = 5;
      out_perm(3) = 6;
      out_perm(4) = 7;
      out_perm(5) = 8;
      out_perm(6) = 9;
      out_perm(7) = 10;
      out_perm(8) = 11;
      out_perm(9) = 12;
      out_perm(10) = 16;
      out_perm(11) = 17;
      out_perm(12) = 18;
      out_perm(13) = 22;
      out_perm(14) = 23;
      out_perm(15) = 24;
      out_perm(16) = 1;
      out_perm(17) = 2;
      out_perm(18) = 3;
      out_perm(19) = 13;
      out_perm(20) = 14;
      out_perm(21) = 15;
      out_perm(22) = 19;
      out_perm(23) = 20;
      out_perm(24) = 21;
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
      out_unperm = zeros(24, 1);
      out_unperm(1) = 16;
      out_unperm(2) = 17;
      out_unperm(3) = 18;
      out_unperm(4) = 1;
      out_unperm(5) = 2;
      out_unperm(6) = 3;
      out_unperm(7) = 4;
      out_unperm(8) = 5;
      out_unperm(9) = 6;
      out_unperm(10) = 7;
      out_unperm(11) = 8;
      out_unperm(12) = 9;
      out_unperm(13) = 19;
      out_unperm(14) = 20;
      out_unperm(15) = 21;
      out_unperm(16) = 10;
      out_unperm(17) = 11;
      out_unperm(18) = 12;
      out_unperm(19) = 22;
      out_unperm(20) = 23;
      out_unperm(21) = 24;
      out_unperm(22) = 13;
      out_unperm(23) = 14;
      out_unperm(24) = 15;
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
end % Test_06

% That's All Folks!
