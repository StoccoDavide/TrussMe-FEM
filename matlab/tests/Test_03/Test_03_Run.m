%% Test_03 Run: 15 x 15 matrix

% Clear workspace
clear all; %#ok<CLALL>
close all;
clc;

% Test title
disp('TrussMe - Test_03: Running...');

% Create an object
obj = Test_03();

% Evaluate stiffness matrices
K_ff = obj.K_ff();
K_sf = obj.K_sf();
K_fs = obj.K_fs();
K_ss = obj.K_ss();
K    = obj.K();

% Evaluate displacements
d_f = obj.d_f();
d_s = obj.d_s();
d   = obj.d();

% Evaluate force vectors
f_f = obj.f_f();
f_s = obj.f_s();
f_r = obj.f_r();
f   = obj.f();

% Compute displacements
d_fc = obj.compute_d_f([], []);
d_c  = obj.compute_d([], []);

% Compute force vectors
f_sc = obj.compute_f_s([], []);
f_c  = obj.compute_f([], []);

% Check size of matrices and vectors
obj.check_size([], []);

% Check symmetry of stiffness matrices
obj.check_symmetry([], []);

% Check deformation solution
assert(norm(d_c - d) < 2.0e-04, 'TrussMe - Test_03: Failed to compute displacements.');

% Check force solution
assert(norm(f_c - f + [f_r; 0; 0; 0; 0]) < 1.0e-00, 'TrussMe - Test_01: Failed to compute forces.');

% Test passed
disp('TrussMe - Test_03: Passed.');

% That's All Folks!
