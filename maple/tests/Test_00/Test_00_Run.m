%% Test_00 Run

% Clear workspace
clear all; %#ok<CLALL>
close all;
clc;

% Test title
disp('TrussMe - Test_00: Running...');

% Create an object
obj = Test_00();

% Evaluate states
x = [];

% Evaluate veiling variables
v = obj.v(x);

% Evaluate stiffness matrices
K_ff = obj.K_ff(x, v);
K_sf = obj.K_sf(x, v);
K_fs = obj.K_fs(x, v);
K_ss = obj.K_ss(x, v);
K    = obj.K(x, v);

% Evaluate displacements
d_f = obj.d_f(x, v);
d_s = obj.d_s(x, v);
d   = obj.d(x, v);

% Evaluate force vectors
f_f = obj.f_f(x, v);
f_s = obj.f_s(x, v);
f_r = obj.f_r(x, v);
f   = obj.f(x, v);

% Compute displacements
d_fc = obj.compute_d_f(x, v);
d_c  = obj.compute_d(x, v);

% Compute force vectors
f_sc = obj.compute_f_s(x, v);
f_c  = obj.compute_f(x, v);

% Check size of matrices and vectors
obj.sanity_check(x, v);

% Check deformation solution
assert(norm(d_c - d) < 1.0e-08, 'TrussMe - Test_00: Failed to compute displacements.');

% Check force solution
assert(norm(f_c - f) < 1.0e-08, 'TrussMe - Test_00: Failed to compute forces.');

% Test passed
disp('TrussMe - Test_00: Passed.');

% That's All Folks!
