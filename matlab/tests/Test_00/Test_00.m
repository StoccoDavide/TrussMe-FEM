% Test_00: 10 x 10 matrix
%
classdef Test_00 < TrussMe.System
  %
  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function this = Test_00()
      this@TrussMe.System([]);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = K( ~, ~, ~ )
      out = 1.0e+06 * [ ...
         0.6667,  0.0000, -0.6667,  0.0000,  0.0000, ...
         0.0000,  0.0000,  0.0000,  0.0000,  0.0000; ...
         0.0000,  0.5000,  0.0000,  0.0000,  0.0000, ...
        -0.5000,  0.0000,  0.0000,  0.0000,  0.0000; ...
        -0.6667,  0.0000,  0.9547,  0.0000, -0.1440, ...
        -0.1920,  0.0000,  0.0000, -0.1440,  0.1920; ...
         0.0000,  0.0000,  0.0000,  1.0120, -0.1920, ...
        -0.2560,  0.0000, -0.5000,  0.1920, -0.2560; ...
         0.0000,  0.0000, -0.1440, -0.1920,  0.8107, ...
         0.1920, -0.6667,  0.0000,  0.0000,  0.0000; ...
         0.0000, -0.5000, -0.1920, -0.2560,  0.1920, ...
         0.7560,  0.0000,  0.0000,  0.0000,  0.0000; ...
         0.0000,  0.0000,  0.0000,  0.0000, -0.6667, ...
         0.0000,  1.3333,  0.0000, -0.6667,  0.0000; ...
         0.0000,  0.0000,  0.0000, -0.5000,  0.0000, ...
         0.0000,  0.0000,  0.5000,  0.0000,  0.0000; ...
         0.0000,  0.0000, -0.1440,  0.1920,  0.0000, ...
         0.0000, -0.6667,  0.0000,  0.8107, -0.1920; ...
         0.0000,  0.0000,  0.1920, -0.2560,  0.0000, ...
         0.0000,  0.0000,  0.0000, -0.1920,  0.2560; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = K_ff( ~, ~, ~ )
      out = 1.0e+06 * [ ...
         0.9547,  0.0000, -0.1920,  0.0000,  0.0000; ...
         0.0000,  1.0120, -0.2560,  0.0000, -0.5000; ...
        -0.1920, -0.2560,  0.7560,  0.0000,  0.0000; ...
         0.0000,  0.0000,  0.0000,  1.3333,  0.0000; ...
         0.0000, -0.5000,  0.0000,  0.0000,  0.5000; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = K_fs( ~, ~, ~ )
      out = 1.0e+06 * [ ...
        -0.6667,  0.0000, -0.1440, -0.1440,  0.1920; ...
         0.0000,  0.0000, -0.1920,  0.1920, -0.2560; ...
         0.0000, -0.5000,  0.1920,  0.0000,  0.0000; ...
         0.0000,  0.0000, -0.6667, -0.6667,  0.0000; ...
         0.0000,  0.0000,  0.0000,  0.0000,  0.0000; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = K_sf( ~, ~, ~ )
      out = 1.0e+06 * [ ...
        -0.6667,  0.0000,  0.0000,  0.0000,  0.0000; ...
         0.0000,  0.0000, -0.5000,  0.0000,  0.0000; ...
        -0.1440, -0.1920,  0.1920, -0.6667,  0.0000; ...
        -0.1440,  0.1920,  0.0000, -0.6667,  0.0000; ...
         0.1920, -0.2560,  0.0000,  0.0000,  0.0000; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = K_ss( ~, ~, ~ )
      out = 1.0e+06 * [ ...
        0.6667,  0.0000,  0.0000,  0.0000,  0.0000; ...
        0.0000,  0.5000,  0.0000,  0.0000,  0.0000; ...
        0.0000,  0.0000,  0.8107,  0.0000,  0.0000; ...
        0.0000,  0.0000,  0.0000,  0.8107, -0.1920; ...
        0.0000,  0.0000,  0.0000, -0.1920,  0.2560; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = d( ~, ~, ~ )
      out = [
          0.0000; ...
          0.0000; ...
          0.0095; ...
         -0.0221; ...
          0.0000; ...
         -0.0051; ...
          0.0000; ...
         -0.0421; ...
          0.0000; ...
          0.0000; ...
       ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = d_f( ~, ~, ~ )
      out = [ ...
         0.0095; ...
        -0.0221; ...
        -0.0051; ...
         0.0000; ...
        -0.0421; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = d_s( ~, ~, ~ )
      out = [ ...
        0.0000; ...
        0.0000; ...
        0.0000; ...
        0.0000; ...
        0.0000; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = f( ~, ~, ~ )
      out = 1.0e+04 * [ ...
        -0.6303; ...
         0.2536; ...
         1.0000; ...
         0.0000; ...
         0.1902; ...
         0.0000; ...
         0.0000; ...
        -1.0000; ...
        -0.5598; ...
         0.7464; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = f_f( ~, ~, ~ )
      out = 1.0e+04 * [ ...
         1.0000; ...
         0.0000; ...
         0.0000; ...
         0.0000; ...
        -1.0000; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = f_s( ~, ~, ~ )
      out = 1.0e+04 * [ ...
        -0.6303
         0.2536
         0.1902
        -0.5598
         0.7464
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = f_r( ~, ~, ~ )
      out = [ ...
        0.0000; ...
        0.0000; ...
        0.0000; ...
        0.0000; ...
        0.0000; ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = perm( ~ )
      out = [3, 4, 6, 7, 8, 1, 2, 5, 9, 10];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = unperm( ~ )
      out = [6, 7, 1, 2, 8, 3, 4, 5, 9, 10];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out = v( ~, ~ )
      out = [];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
  %
end

% That's All Folks!