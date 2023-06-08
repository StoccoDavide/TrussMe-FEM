% +--------------------------------------------------------------------------+
% | 'TrussMe' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023 |
% | Current version authors:                                                 |
% |   Davide Stocco and Matteo Larcher.                                      |
% +--------------------------------------------------------------------------+

% Matlab generated code for system: Test_00
% This file has been automatically generated by TrussMe.
% DISCLAIMER: If you need to edit it, do it wisely!

classdef Test_00 < TrussMe.System
  %
  % Test_00 class
  %
  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function this = Test_00( varargin )
      % Class constructor.

      % User data
      if (nargin == 0)
        data.N1_y = 4;
      elseif (nargin == 1 && isstruct(varargin{1}))
        data = varargin{1};
      elseif (nargin == 1)
        data.N1_y = varargin{1};
      else
        error('wrong number of input arguments.');
      end

      % Call superclass constructor
      this = this@TrussMe.System(data);
    end % Test_00
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K = K( this, ~, ~ )
      % Evaluate the stiffness matrix K.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      t1 = 4 - N1_y;
      t2 = t1 ^ 2;
      t3 = 9 + t2;
      t4 = sqrt(t3);
      t6 = 0.1e1 / t4 / t3;
      out_1_1 = 0.1800e8 * t6;
      out_2_1 = 0.600e7 * t1 * t6;
      out_7_1 = -out_1_1;
      out_8_1 = -out_2_1;
      out_1_2 = out_2_1;
      t9 = 0.200e7 * t2 * t6;
      t10 = N1_y ^ 2;
      t11 = sqrt(t10);
      t13 = 0.200e7 / t11;
      out_2_2 = t9 + t13;
      out_7_2 = out_8_1;
      out_8_2 = -t9;
      out_14_2 = -t13;
      out_6_6 = 1;
      out_1_7 = out_7_1;
      out_2_7 = out_7_2;
      out_7_7 = out_1_1 + 0.288000e6;
      out_8_7 = out_1_2;
      out_13_7 = -0.144000e6;
      out_14_7 = -0.192000e6;
      out_25_7 = -0.144000e6;
      out_26_7 = 0.192000e6;
      out_1_8 = out_2_7;
      out_2_8 = out_8_2;
      out_7_8 = out_8_7;
      out_8_8 = t9 + 0.1012000e7;
      out_13_8 = -0.192000e6;
      out_14_8 = -0.256000e6;
      out_20_8 = -0.500000e6;
      out_25_8 = 0.192000e6;
      out_26_8 = -0.256000e6;
      out_9_9 = 1;
      out_10_10 = 1;
      out_11_11 = 1;
      out_12_12 = 1;
      out_7_13 = -0.144000e6;
      out_8_13 = -0.192000e6;
      out_13_13 = 0.810666666699999943e6;
      out_14_13 = 0.192000e6;
      out_19_13 = -0.666666666699999943e6;
      out_2_14 = out_14_2;
      out_7_14 = -0.192000e6;
      out_8_14 = -0.256000e6;
      out_13_14 = 0.192000e6;
      out_14_14 = t13 + 0.256000e6;
      out_18_18 = 1;
      out_13_19 = -0.666666666699999943e6;
      out_19_19 = 0.133333333339999989e7;
      out_25_19 = -0.666666666699999943e6;
      out_8_20 = -0.500000e6;
      out_20_20 = 0.500000e6;
      out_21_21 = 1;
      out_22_22 = 1;
      out_23_23 = 1;
      out_24_24 = 1;
      out_7_25 = -0.144000e6;
      out_8_25 = 0.192000e6;
      out_19_25 = -0.666666666699999943e6;
      out_25_25 = 0.810666666699999943e6;
      out_26_25 = -0.192000e6;
      out_7_26 = 0.192000e6;
      out_8_26 = -0.256000e6;
      out_25_26 = -0.192000e6;
      out_26_26 = 0.256000e6;
      out_30_30 = 1;

      % Store outputs
      out_K = zeros(30, 30);
      out_K(1, 1) = out_1_1;
      out_K(2, 1) = out_2_1;
      out_K(7, 1) = out_7_1;
      out_K(8, 1) = out_8_1;
      out_K(1, 2) = out_1_2;
      out_K(2, 2) = out_2_2;
      out_K(7, 2) = out_7_2;
      out_K(8, 2) = out_8_2;
      out_K(14, 2) = out_14_2;
      out_K(6, 6) = out_6_6;
      out_K(1, 7) = out_1_7;
      out_K(2, 7) = out_2_7;
      out_K(7, 7) = out_7_7;
      out_K(8, 7) = out_8_7;
      out_K(13, 7) = out_13_7;
      out_K(14, 7) = out_14_7;
      out_K(25, 7) = out_25_7;
      out_K(26, 7) = out_26_7;
      out_K(1, 8) = out_1_8;
      out_K(2, 8) = out_2_8;
      out_K(7, 8) = out_7_8;
      out_K(8, 8) = out_8_8;
      out_K(13, 8) = out_13_8;
      out_K(14, 8) = out_14_8;
      out_K(20, 8) = out_20_8;
      out_K(25, 8) = out_25_8;
      out_K(26, 8) = out_26_8;
      out_K(9, 9) = out_9_9;
      out_K(10, 10) = out_10_10;
      out_K(11, 11) = out_11_11;
      out_K(12, 12) = out_12_12;
      out_K(7, 13) = out_7_13;
      out_K(8, 13) = out_8_13;
      out_K(13, 13) = out_13_13;
      out_K(14, 13) = out_14_13;
      out_K(19, 13) = out_19_13;
      out_K(2, 14) = out_2_14;
      out_K(7, 14) = out_7_14;
      out_K(8, 14) = out_8_14;
      out_K(13, 14) = out_13_14;
      out_K(14, 14) = out_14_14;
      out_K(18, 18) = out_18_18;
      out_K(13, 19) = out_13_19;
      out_K(19, 19) = out_19_19;
      out_K(25, 19) = out_25_19;
      out_K(8, 20) = out_8_20;
      out_K(20, 20) = out_20_20;
      out_K(21, 21) = out_21_21;
      out_K(22, 22) = out_22_22;
      out_K(23, 23) = out_23_23;
      out_K(24, 24) = out_24_24;
      out_K(7, 25) = out_7_25;
      out_K(8, 25) = out_8_25;
      out_K(19, 25) = out_19_25;
      out_K(25, 25) = out_25_25;
      out_K(26, 25) = out_26_25;
      out_K(7, 26) = out_7_26;
      out_K(8, 26) = out_8_26;
      out_K(25, 26) = out_25_26;
      out_K(26, 26) = out_26_26;
      out_K(30, 30) = out_30_30;
    end % K
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_ff = K_ff( this, ~, ~ )
      % Evaluate the stiffness matrix K_ff.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      out_1_1 = 1;
      t1 = 4 - N1_y;
      t2 = t1 ^ 2;
      t3 = 9 + t2;
      t4 = sqrt(t3);
      t6 = 0.1e1 / t4 / t3;
      out_2_2 = 0.1800e8 * t6 + 0.288000e6;
      out_3_2 = 0.600e7 * t1 * t6;
      out_8_2 = -0.192000e6;
      out_2_3 = out_3_2;
      out_3_3 = 0.200e7 * t2 * t6 + 0.1012000e7;
      out_8_3 = -0.256000e6;
      out_11_3 = -0.500000e6;
      out_4_4 = 1;
      out_5_5 = 1;
      out_6_6 = 1;
      out_7_7 = 1;
      out_2_8 = -0.192000e6;
      out_3_8 = -0.256000e6;
      t11 = N1_y ^ 2;
      t12 = sqrt(t11);
      out_8_8 = 0.200e7 / t12 + 0.256000e6;
      out_9_9 = 1;
      out_10_10 = 0.133333333339999989e7;
      out_3_11 = -0.500000e6;
      out_11_11 = 0.500000e6;
      out_12_12 = 1;
      out_13_13 = 1;
      out_14_14 = 1;
      out_15_15 = 1;
      out_16_16 = 1;

      % Store outputs
      out_K_ff = zeros(16, 16);
      out_K_ff(1, 1) = out_1_1;
      out_K_ff(2, 2) = out_2_2;
      out_K_ff(3, 2) = out_3_2;
      out_K_ff(8, 2) = out_8_2;
      out_K_ff(2, 3) = out_2_3;
      out_K_ff(3, 3) = out_3_3;
      out_K_ff(8, 3) = out_8_3;
      out_K_ff(11, 3) = out_11_3;
      out_K_ff(4, 4) = out_4_4;
      out_K_ff(5, 5) = out_5_5;
      out_K_ff(6, 6) = out_6_6;
      out_K_ff(7, 7) = out_7_7;
      out_K_ff(2, 8) = out_2_8;
      out_K_ff(3, 8) = out_3_8;
      out_K_ff(8, 8) = out_8_8;
      out_K_ff(9, 9) = out_9_9;
      out_K_ff(10, 10) = out_10_10;
      out_K_ff(3, 11) = out_3_11;
      out_K_ff(11, 11) = out_11_11;
      out_K_ff(12, 12) = out_12_12;
      out_K_ff(13, 13) = out_13_13;
      out_K_ff(14, 14) = out_14_14;
      out_K_ff(15, 15) = out_15_15;
      out_K_ff(16, 16) = out_16_16;
    end % K_ff
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_fs = K_fs( this, ~, ~ )
      % Evaluate the stiffness matrix K_fs.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      t1 = 4 - N1_y;
      t2 = t1 ^ 2;
      t3 = 9 + t2;
      t4 = sqrt(t3);
      t6 = 0.1e1 / t4 / t3;
      out_2_1 = -0.1800e8 * t6;
      out_3_1 = -0.600e7 * t1 * t6;
      out_2_2 = out_3_1;
      out_3_2 = -0.200e7 * t2 * t6;
      t12 = N1_y ^ 2;
      t13 = sqrt(t12);
      out_8_2 = -0.200e7 / t13;
      out_2_6 = -0.144000e6;
      out_3_6 = -0.192000e6;
      out_8_6 = 0.192000e6;
      out_10_6 = -0.666666666699999943e6;
      out_2_10 = -0.144000e6;
      out_3_10 = 0.192000e6;
      out_10_10 = -0.666666666699999943e6;
      out_2_11 = 0.192000e6;
      out_3_11 = -0.256000e6;

      % Store outputs
      out_K_fs = zeros(16, 14);
      out_K_fs(2, 1) = out_2_1;
      out_K_fs(3, 1) = out_3_1;
      out_K_fs(2, 2) = out_2_2;
      out_K_fs(3, 2) = out_3_2;
      out_K_fs(8, 2) = out_8_2;
      out_K_fs(2, 6) = out_2_6;
      out_K_fs(3, 6) = out_3_6;
      out_K_fs(8, 6) = out_8_6;
      out_K_fs(10, 6) = out_10_6;
      out_K_fs(2, 10) = out_2_10;
      out_K_fs(3, 10) = out_3_10;
      out_K_fs(10, 10) = out_10_10;
      out_K_fs(2, 11) = out_2_11;
      out_K_fs(3, 11) = out_3_11;
    end % K_fs
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_sf = K_sf( this, ~, ~ )
      % Evaluate the stiffness matrix K_sf.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      t1 = 4 - N1_y;
      t2 = t1 ^ 2;
      t3 = 9 + t2;
      t4 = sqrt(t3);
      t6 = 0.1e1 / t4 / t3;
      out_1_2 = -0.1800e8 * t6;
      out_2_2 = -0.600e7 * t1 * t6;
      out_6_2 = -0.144000e6;
      out_10_2 = -0.144000e6;
      out_11_2 = 0.192000e6;
      out_1_3 = out_2_2;
      out_2_3 = -0.200e7 * t2 * t6;
      out_6_3 = -0.192000e6;
      out_10_3 = 0.192000e6;
      out_11_3 = -0.256000e6;
      t12 = N1_y ^ 2;
      t13 = sqrt(t12);
      out_2_8 = -0.200e7 / t13;
      out_6_8 = 0.192000e6;
      out_6_10 = -0.666666666699999943e6;
      out_10_10 = -0.666666666699999943e6;

      % Store outputs
      out_K_sf = zeros(14, 16);
      out_K_sf(1, 2) = out_1_2;
      out_K_sf(2, 2) = out_2_2;
      out_K_sf(6, 2) = out_6_2;
      out_K_sf(10, 2) = out_10_2;
      out_K_sf(11, 2) = out_11_2;
      out_K_sf(1, 3) = out_1_3;
      out_K_sf(2, 3) = out_2_3;
      out_K_sf(6, 3) = out_6_3;
      out_K_sf(10, 3) = out_10_3;
      out_K_sf(11, 3) = out_11_3;
      out_K_sf(2, 8) = out_2_8;
      out_K_sf(6, 8) = out_6_8;
      out_K_sf(6, 10) = out_6_10;
      out_K_sf(10, 10) = out_10_10;
    end % K_sf
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_K_ss = K_ss( this, ~, ~ )
      % Evaluate the stiffness matrix K_ss.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      t1 = 4 - N1_y;
      t2 = t1 ^ 2;
      t3 = 9 + t2;
      t4 = sqrt(t3);
      t6 = 0.1e1 / t4 / t3;
      out_1_1 = 0.1800e8 * t6;
      out_2_1 = 0.600e7 * t1 * t6;
      out_1_2 = out_2_1;
      t10 = N1_y ^ 2;
      t11 = sqrt(t10);
      out_2_2 = 0.200e7 * t2 * t6 + 0.200e7 / t11;
      out_6_6 = 0.810666666699999943e6;
      out_10_10 = 0.810666666699999943e6;
      out_11_10 = -0.192000e6;
      out_10_11 = -0.192000e6;
      out_11_11 = 0.256000e6;

      % Store outputs
      out_K_ss = zeros(14, 14);
      out_K_ss(1, 1) = out_1_1;
      out_K_ss(2, 1) = out_2_1;
      out_K_ss(1, 2) = out_1_2;
      out_K_ss(2, 2) = out_2_2;
      out_K_ss(6, 6) = out_6_6;
      out_K_ss(10, 10) = out_10_10;
      out_K_ss(11, 10) = out_11_10;
      out_K_ss(10, 11) = out_10_11;
      out_K_ss(11, 11) = out_11_11;
    end % K_ss
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_d = d( this, ~, in_2 )
      % Evaluate the deformation vector d.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      V_1 = in_2(1);
      V_2 = in_2(2);
      V_3 = in_2(3);
      V_4 = in_2(4);
      V_5 = in_2(5);
      V_6 = in_2(6);
      V_7 = in_2(7);
      V_8 = in_2(8);
      V_9 = in_2(9);
      V_10 = in_2(10);

      % Evaluate function
      out_7 = V_10;
      out_8 = -0.1e1 / (-0.600e7 * V_3 - 0.750000000000000000e0 * V_4 + 0.375000e6 + V_6 * V_1 + V_2 * V_6 * 0.4500000e7) * V_7;
      out_14 = V_9;
      out_20 = -V_8;

      % Store outputs
      out_d = zeros(30, 1);
      out_d(7) = out_7;
      out_d(8) = out_8;
      out_d(14) = out_14;
      out_d(20) = out_20;
    end % d
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_d_f = d_f( this, ~, in_2 )
      % Evaluate the deformation vector d_f.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      V_1 = in_2(1);
      V_2 = in_2(2);
      V_3 = in_2(3);
      V_4 = in_2(4);
      V_5 = in_2(5);
      V_6 = in_2(6);
      V_7 = in_2(7);
      V_8 = in_2(8);
      V_9 = in_2(9);
      V_10 = in_2(10);

      % Evaluate function
      out_2 = V_10;
      out_3 = -0.1e1 / (-0.600e7 * V_3 - 0.750000000000000000e0 * V_4 + 0.375000e6 + V_6 * V_1 + V_2 * V_6 * 0.4500000e7) * V_7;
      out_8 = V_9;
      out_11 = -V_8;

      % Store outputs
      out_d_f = zeros(16, 1);
      out_d_f(2) = out_2;
      out_d_f(3) = out_3;
      out_d_f(8) = out_8;
      out_d_f(11) = out_11;
    end % d_f
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_d_s = d_s( this, ~, ~ )
      % Evaluate the deformation vector d_s.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_d_s = zeros(14, 1);
    end % d_s
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f = f( this, ~, in_2 )
      % Evaluate the force vector f.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      V_1 = in_2(1);
      V_2 = in_2(2);
      V_3 = in_2(3);
      V_4 = in_2(4);
      V_5 = in_2(5);
      V_6 = in_2(6);
      V_7 = in_2(7);
      V_8 = in_2(8);
      V_9 = in_2(9);
      V_10 = in_2(10);

      % Evaluate function
      t1 = 4 - N1_y;
      t2 = t1 ^ 2;
      t3 = 9 + t2;
      t4 = sqrt(t3);
      t6 = 0.1e1 / t4 / t3;
      t9 = t1 * t6;
      t17 = 0.1e1 / (-0.600e7 * V_3 - 0.750000000000000000e0 * V_4 + 0.375000e6 + V_6 * V_1 + V_2 * V_6 * 0.4500000e7) * V_7;
      out_1 = -0.1800e8 * V_10 * t6 + 0.600e7 * t17 * t9;
      t25 = N1_y ^ 2;
      t26 = sqrt(t25);
      out_2 = -0.600e7 * V_10 * t9 + 0.200e7 * t17 * t2 * t6 - 0.200e7 * V_9 / t26;
      out_7 = 10000;
      t30 = 0.144000e6 * V_10;
      out_13 = -t30 + t17 * 0.192000e6 + V_9 * 0.192000e6;
      out_20 = -10000;
      out_25 = -t30 - 0.192000e6 * t17;
      out_26 = V_10 * 0.192000e6 + t17 * 0.256000e6;

      % Store outputs
      out_f = zeros(30, 1);
      out_f(1) = out_1;
      out_f(2) = out_2;
      out_f(7) = out_7;
      out_f(13) = out_13;
      out_f(20) = out_20;
      out_f(25) = out_25;
      out_f(26) = out_26;
    end % f
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f_f = f_f( this, ~, ~ )
      % Evaluate the force vector f_f.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      out_2 = 10000;
      out_11 = -10000;

      % Store outputs
      out_f_f = zeros(16, 1);
      out_f_f(2) = out_2;
      out_f_f(11) = out_11;
    end % f_f
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f_s = f_s( this, ~, in_2 )
      % Evaluate the force vector f_s.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      V_1 = in_2(1);
      V_2 = in_2(2);
      V_3 = in_2(3);
      V_4 = in_2(4);
      V_5 = in_2(5);
      V_6 = in_2(6);
      V_7 = in_2(7);
      V_8 = in_2(8);
      V_9 = in_2(9);
      V_10 = in_2(10);

      % Evaluate function
      t1 = 4 - N1_y;
      t2 = t1 ^ 2;
      t3 = 9 + t2;
      t4 = sqrt(t3);
      t6 = 0.1e1 / t4 / t3;
      t9 = t1 * t6;
      t17 = 0.1e1 / (-0.600e7 * V_3 - 0.750000000000000000e0 * V_4 + 0.375000e6 + V_6 * V_1 + V_2 * V_6 * 0.4500000e7) * V_7;
      out_1 = -0.1800e8 * V_10 * t6 + 0.600e7 * t17 * t9;
      t25 = N1_y ^ 2;
      t26 = sqrt(t25);
      out_2 = -0.600e7 * V_10 * t9 + 0.200e7 * t17 * t2 * t6 - 0.200e7 * V_9 / t26;
      t30 = 0.144000e6 * V_10;
      out_6 = -t30 + t17 * 0.192000e6 + V_9 * 0.192000e6;
      out_10 = -t30 - 0.192000e6 * t17;
      out_11 = V_10 * 0.192000e6 + t17 * 0.256000e6;

      % Store outputs
      out_f_s = zeros(14, 1);
      out_f_s(1) = out_1;
      out_f_s(2) = out_2;
      out_f_s(6) = out_6;
      out_f_s(10) = out_10;
      out_f_s(11) = out_11;
    end % f_s
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_f_r = f_r( this, ~, ~ )
      % Evaluate the force vector f_r.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      % No body

      % Store outputs
      out_f_r = zeros(14, 1);
    end % f_r
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_perm = perm( this )
      % Evaluate the permutation vector.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      out_1 = 6;
      out_2 = 7;
      out_3 = 8;
      out_4 = 9;
      out_5 = 10;
      out_6 = 11;
      out_7 = 12;
      out_8 = 14;
      out_9 = 18;
      out_10 = 19;
      out_11 = 20;
      out_12 = 21;
      out_13 = 22;
      out_14 = 23;
      out_15 = 24;
      out_16 = 30;
      out_17 = 1;
      out_18 = 2;
      out_19 = 3;
      out_20 = 4;
      out_21 = 5;
      out_22 = 13;
      out_23 = 15;
      out_24 = 16;
      out_25 = 17;
      out_26 = 25;
      out_27 = 26;
      out_28 = 27;
      out_29 = 28;
      out_30 = 29;

      % Store outputs
      out_perm = zeros(30, 1);
      out_perm(1) = out_1;
      out_perm(2) = out_2;
      out_perm(3) = out_3;
      out_perm(4) = out_4;
      out_perm(5) = out_5;
      out_perm(6) = out_6;
      out_perm(7) = out_7;
      out_perm(8) = out_8;
      out_perm(9) = out_9;
      out_perm(10) = out_10;
      out_perm(11) = out_11;
      out_perm(12) = out_12;
      out_perm(13) = out_13;
      out_perm(14) = out_14;
      out_perm(15) = out_15;
      out_perm(16) = out_16;
      out_perm(17) = out_17;
      out_perm(18) = out_18;
      out_perm(19) = out_19;
      out_perm(20) = out_20;
      out_perm(21) = out_21;
      out_perm(22) = out_22;
      out_perm(23) = out_23;
      out_perm(24) = out_24;
      out_perm(25) = out_25;
      out_perm(26) = out_26;
      out_perm(27) = out_27;
      out_perm(28) = out_28;
      out_perm(29) = out_29;
      out_perm(30) = out_30;
    end % perm
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_unperm = unperm( this )
      % Evaluate the unpermutation vector.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      out_1 = 17;
      out_2 = 18;
      out_3 = 19;
      out_4 = 20;
      out_5 = 21;
      out_6 = 1;
      out_7 = 2;
      out_8 = 3;
      out_9 = 4;
      out_10 = 5;
      out_11 = 6;
      out_12 = 7;
      out_13 = 22;
      out_14 = 8;
      out_15 = 23;
      out_16 = 24;
      out_17 = 25;
      out_18 = 9;
      out_19 = 10;
      out_20 = 11;
      out_21 = 12;
      out_22 = 13;
      out_23 = 14;
      out_24 = 15;
      out_25 = 26;
      out_26 = 27;
      out_27 = 28;
      out_28 = 29;
      out_29 = 30;
      out_30 = 16;

      % Store outputs
      out_unperm = zeros(30, 1);
      out_unperm(1) = out_1;
      out_unperm(2) = out_2;
      out_unperm(3) = out_3;
      out_unperm(4) = out_4;
      out_unperm(5) = out_5;
      out_unperm(6) = out_6;
      out_unperm(7) = out_7;
      out_unperm(8) = out_8;
      out_unperm(9) = out_9;
      out_unperm(10) = out_10;
      out_unperm(11) = out_11;
      out_unperm(12) = out_12;
      out_unperm(13) = out_13;
      out_unperm(14) = out_14;
      out_unperm(15) = out_15;
      out_unperm(16) = out_16;
      out_unperm(17) = out_17;
      out_unperm(18) = out_18;
      out_unperm(19) = out_19;
      out_unperm(20) = out_20;
      out_unperm(21) = out_21;
      out_unperm(22) = out_22;
      out_unperm(23) = out_23;
      out_unperm(24) = out_24;
      out_unperm(25) = out_25;
      out_unperm(26) = out_26;
      out_unperm(27) = out_27;
      out_unperm(28) = out_28;
      out_unperm(29) = out_29;
      out_unperm(30) = out_30;
    end % unperm
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out_v = v( this, ~ )
      % Evaluate the veiling vector.

      % Extract properties
      N1_y = this.m_data.N1_y;

      % Extract inputs
      % No inputs

      % Evaluate function
      t1 = N1_y ^ 2;
      t3 = t1 - 8 * N1_y + 25;
      t4 = sqrt(t3);
      t5 = t4 * t3;
      t8 = 0.1e1 / t5;
      V_1 = t8 * (0.1800e8 + t5 * 0.288000e6);
      V_2 = 0.1000000000e1 * t8 * (-4 + N1_y);
      V_3 = V_2;
      V_4 = t8 * (0.200e7 * t1 - 0.1600e8 * N1_y + 0.3200e8 + t5 * 0.1012000e7);
      t15 = sqrt(t1);
      V_5 = 0.1e1 / t15 * (0.200e7 + t15 * 0.256000e6);
      t19 = V_5 * V_4;
      t23 = V_5 * V_2;
      t26 = 0.1e1 / (0.192000e6 + 0.2343750000e2 * t23);
      V_6 = t26 * (-0.256000e6 + 0.3906250000e-5 * t19 - 0.1953125000e1 * V_5);
      V_7 = t26 * (0.3906250000e-1 * V_5 * V_1 - 0.234375e6 * t23 - 0.3360000000e10);
      t33 = V_6 * V_1;
      t35 = V_2 * V_6;
      t42 = 0.1e1 / (-0.6000000e7 * V_3 - 0.7500000000e0 * V_4 + 0.375000e6 + t33 + 0.4500000e7 * t35);
      V_8 = t42 * (V_7 - 0.120000e6 * V_3 - 0.1500000000e-1 * V_4 + 0.7500e4 + 0.2000000000e-1 * t33 + 0.90000e5 * t35);
      t43 = V_7 * V_4;
      t54 = V_2 ^ 2;
      t56 = V_6 * V_5;
      t63 = -0.7500000000e0 * t43 - 0.9155273438e-4 * t23 * t43 + 0.375000e6 * V_7 + 0.4577636719e2 * V_7 * V_5 * V_2 + 0.4500000e7 * V_7 * V_2 * V_6 + 0.5493164062e3 * t56 * t54 * V_7 - 0.4500000000e11 * V_3 - 0.5625e4 * V_4 + 0.2812500000e10 + 0.7500e4 * t33 + 0.3375000000e11 * t35;
      V_9 = t42 * t26 * t63;
      t65 = V_7 * V_6;
      V_10 = t42 * t26 * (-0.192000e6 * t65 - 0.2343750000e2 * t23 * t65 - 0.234375e6 * V_5 * V_3 - 0.2929687500e-1 * t19 + 0.1464843750e5 * V_5 + 0.3906250000e-1 * V_1 * t56 + 0.1757812500e6 * V_6 * t23);

      % Store outputs
      out_v = zeros(10, 1);
      out_v(1) = V_1;
      out_v(2) = V_2;
      out_v(3) = V_3;
      out_v(4) = V_4;
      out_v(5) = V_5;
      out_v(6) = V_6;
      out_v(7) = V_7;
      out_v(8) = V_8;
      out_v(9) = V_9;
      out_v(10) = V_10;
    end % v
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
end % Test_00

% That's All Folks!