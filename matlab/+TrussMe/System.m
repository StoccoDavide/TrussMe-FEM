%
%> Class container for the system
%
classdef System < handle
  %
  properties (SetAccess = protected, Hidden = true)
    %
    %> System data.
    %
    m_data;
    %
  end
  %
  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Class constructor for a system.:
    %>
    %> \param t_data The system data.
    %>
    %> \return The system.
    %
    function this = System( t_data )
      this.m_data = t_data;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Get the system data.
    %>
    %> \return The system data.
    %
    function out = get_data( this )
      out = this.m_data;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Set the system data.
    %>
    %> \param t_data The system data.
    %
    function set_data( this, t_data )
      this.m_data = t_data;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the system stiffness matrix \f$ \mathbf{K} \f$ as:
    %>
    %> \f[
    %>   \mathbf{K} = \left[
    %>     \begin{array}{cc}
    %>       \mathbf{K}_{ff} & \mathbf{K}_{fs} \\
    %>       \mathbf{K}_{sf} & \mathbf{K}_{ss}
    %>     \end{array}
    %>   \right]
    %> \f]
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system stiffness matrix \f$ \mathbf{K} \f$.
    %
    function out = compute_K( this, x, v )
      out = [ ...
        this.K_ff(x, v), this.K_fs(x, v); ...
        this.K_sf(x, v), this.K_ss(x, v); ...
      ];
      out = out(this.unperm(), this.unperm());
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the system deformation vector of free dofs \f$ \mathbf{d}_{f} \f$
    %> as:
    %>
    %> \f[
    %>   \mathbf{d}_{f} = \mathbf{K}_{ff}^{-1} \left(
    %>     \mathbf{f}_{f} - \mathbf{K}_{fs} \mathbf{d}_{s}
    %>   \right)
    %> \f]
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{d}_{f} \f$.
    %
    function out = compute_d_f( this, x, v )
      out = this.K_ff(x, v)\(this.f_f(x, v)-this.K_fs(x, v)*this.d_s(x, v));
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the system deformation vector \f$ \mathbf{d} \f$ as:
    %>
    %> \f[
    %>   \mathbf{d} = \left[
    %>     \begin{array}{c} \mathbf{d}_{f} \\ \mathbf{d}_{s} \end{array}
    %>   \right]
    %> \f]
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{d} \f$.
    %
    function out = compute_d( this, x, v )
      d_s = this.d_s(x, v);
      out = [ ...
        this.K_ff(x, v)\(this.f_f(x, v)-this.K_fs(x, v)*d_s); ...
        d_s ...
      ];
      out = out(this.unperm());
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the system force vector of specified dofs \f$ \mathbf{f}_{s} \f$
    %> as:
    %>
    %> \f[
    %>   \mathbf{f}_{s} = \mathbf{K}_{sf} \mathbf{d}_{f} +
    %>     \mathbf{K}_{ss} \mathbf{d}_{s} - \mathbf{f}_{r}
    %> \f]
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{f}_{s} \f$.
    %
    function out = compute_f_s( this, x, v )
      out = this.K_sf(x, v)*this.compute_d_f(x, v) + ...
            this.K_ss(x, v)*this.d_s(x, v) - this.f_r(x, v);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the system force vector \f$ \mathbf{f} \f$ as:
    %>
    %> \f[
    %>   \mathbf{f} = \left[
    %>     \begin{array}{c} \mathbf{f}_{f} \\ \mathbf{f}_{s} \end{array}
    %>   \right]
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{f} \f$.
    %
    function out = compute_f( this, x, v )
      out = [ ...
        this.f_f(x, v); ...
        this.K_sf(x, v)*this.compute_d_f(x, v) + ...
        this.K_ss(x, v)*this.d_s(x, v) - this.f_r(x, v); ...
      ];
      out = out(this.unperm());
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Internal matrices size check.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return An error is thrown if the sizes are not correct.
    %
    function check_size( this, x, v )

      CMD = 'TrussMe.System.check_size(...): ';

      % Evaluate stiffness matrices
      K_ff = this.K_ff(x, v);
      K_sf = this.K_sf(x, v);
      K_fs = this.K_fs(x, v);
      K_ss = this.K_ss(x, v);
      K    = this.K(x, v);
      K_c  = this.compute_K(x, v);

      % Evaluate displacements
      d_f = this.d_f(x, v);
      d_s = this.d_s(x, v);
      d   = this.d(x, v);

      % Evaluate force vectors
      f_f = this.f_f(x, v);
      f_s = this.f_s(x, v);
      f_r = this.f_r(x, v);
      f   = this.f(x, v);

      % Compute displacements
      d_fc = this.compute_d_f(x, v);
      d_c = this.compute_d(x, v);

      % Compute force vectors
      f_sc = this.compute_f_s(x, v);
      f_c  = this.compute_f(x, v);

      % Check sizes
      assert(size(K_ff, 1) == size(K_ff, 2), [CMD, 'K_ff is not square.']);
      assert(size(K_fs, 1) == size(K_fs, 2), [CMD, 'K_fs is not square.']);
      assert(size(K_sf, 1) == size(K_sf, 2), [CMD, 'K_sf is not square.']);
      assert(size(K_ss, 1) == size(K_ss, 2), [CMD, 'K_ss is not square.']);
      assert(size(K, 1)    == size(K, 2),    [CMD, 'K is not square.']);
      assert(size(K_c, 1)  == size(K_c, 2),  [CMD, 'computed K is not square.']);

      assert(size(K_ff, 1) + size(K_fs, 1) == size(K, 1) && ...
             size(K_ff, 2) + size(K_sf, 2) == size(K, 2), ...
        [CMD, 'K_ff, K_sf and K are not compatible.']);
      assert(size(K_ss, 1) + size(K_sf, 1) == size(K, 1) && ...
             size(K_ss, 2) + size(K_fs, 2) == size(K, 2), ...
        [CMD, 'K_ss, K_fs and K are not compatible.']);
      assert(size(K_c, 1) == size(K, 1) && size(K_c, 2) == size(K, 2), ...
        [CMD, 'computed K and K are not compatible.']);

      assert(size(d_f, 2) == 1 && size(d_f, 1) == size(K_ff, 1), ...
        [CMD, 'd_f and K_ff are not compatible.']);
      assert(size(d_fc, 2) == 1 && size(d_fc, 1) == size(K_ff, 1), ...
        [CMD 'computed d_f and K_ff are not compatible.']);
      assert(size(d_s, 2) == 1 && size(d_s, 1) == size(K_ss, 1), ...
        [CMD 'd_s and K_ss are not compatible.']);
      assert(size(d, 2) == 1 && size(d, 1) == size(K, 1), ...
        [CMD, 'd and K are not compatible.']);
      assert(size(d_c, 2) == 1 && size(d_c, 1) == size(K, 1), ...
        [CMD, 'computed d and K are not compatible.']);

      assert(size(f_f, 2) == 1 && size(f_f, 1) == size(K_ff, 1), ...
        [CMD, 'f_f and K_ff are not compatible.']);
      assert(size(f_s, 2) == 1 && size(f_s, 1) == size(K_ss, 1), ...
        [CMD, 'f_s and K_ss are not compatible.']);
      assert(size(f_sc, 2) == 1 && size(f_sc, 1) == size(K_ss, 1), ...
        [CMD, 'computed f_s and K_ss are not compatible.']);
      assert(size(f_r, 2) == 1 && size(f_r, 1) == size(K_ss, 1), ...
        [CMD, 'f_r and K_ss are not compatible.']);
      assert(size(f, 2) == 1 && size(f, 1) == size(K, 1), ...
        [CMD, 'f and K are not compatible.']);
      assert(size(f_c, 2) == 1 && size(f_c, 1) == size(K, 1), ...
        [CMD, 'computed f and K are not compatible.']);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Check symmetry of stiffness matrices.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return An error is thrown if the matrices are not symmetric.
    %
    function check_symmetry( this, x, v )

      CMD = 'TrussMe.System.check_symmetry(...): ';

      % Evaluate stiffness matrices
      K_ff = this.K_ff(x, v);
      K_ss = this.K_ss(x, v);
      K    = this.K(x, v);
      K_c  = this.compute_K(x, v);

      % Check symmetry
      assert(isequal(K_ff, K_ff'), [CMD, 'K_ff is not symmetric.']);
      assert(isequal(K_ss, K_ss'), [CMD, 'K_ss is not symmetric.']);
      assert(isequal(K, K'),       [CMD, 'K is not symmetric.']);
      assert(isequal(K_c, K_c'),   [CMD, 'computed K is not symmetric.']);
    end
  end
  %
  methods (Abstract)
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system stiffness matrix \f$ \mathbf{K} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system stiffness matrix \f$ \mathbf{K} \f$.
    %
    K( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system stiffness matrix of free-free dofs
    %> \f$ \mathbf{K}_{ff} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system stiffness matrix \f$ \mathbf{K}_{ff} \f$.
    %
    K_ff( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system stiffness matrix of free-specified dofs
    %> \f$ \mathbf{K}_{fs} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system stiffness matrix \f$ \mathbf{K}_{fs} \f$.
    %
    K_fs( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system stiffness matrix of specified-free dofs
    %> \f$ \mathbf{K}_{sf} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system stiffness matrix \f$ \mathbf{K}_{sf} \f$.
    %
    K_sf( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system stiffness matrix of specified-specified dofs
    %> \f$ \mathbf{K}_{ss} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system stiffness matrix \f$ \mathbf{K}_{ss} \f$.
    %
    K_ss( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system deformation vector \f$ \mathbf{d} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{d} \f$.
    %
    d( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system deformation vector of free dofs
    %> \f$ \mathbf{d}_{f} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{d}_{f} \f$.
    %
    d_f( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system deformation vector of specified dofs
    %> \f$ \mathbf{d}_{s} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{d}_{s} \f$.
    %
    d_s( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system force vector \f$ \mathbf{f} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{f} \f$.
    %
    f( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system force vector of free dofs \f$ \mathbf{f}_{f} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{f}_{f} \f$.
    %
    f_f( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system force vector of specified dofs \f$ \mathbf{f}_{s} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{f}_{s} \f$.
    %
    f_s( this, x, v )
    %
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the system remainder force vector \f$ \mathbf{f}_{r} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system remainder force vector \f$ \mathbf{f}_{r} \f$.
    %
    f_r( this, x, v )
    %
    %
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Get the system permutation vector.
    %>
    %> \return The permutation vector.
    %
    perm( this )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Get the unpermutation vector.
    %>
    %> \return The unpermutation vector.
    %
    unperm( this )
    %
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the veils \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %>
    %> \return The Veils \f$ \mathbf{v} \f$..
    %
    v( this, x )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
  %
end

% That's All Folks!
