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
    %> Number of free dofs of the system.
    %
    m_num_free_dofs;
    %
    %> Number of specified dofs of the system.
    %
    m_num_spec_dofs;
    %
    %> Number of veiling variables \mathbf{v} of the system.
    %
    m_num_veils;
    %
    %
    %> Number of parameters \mathbf{x} of the system.
    %
    m_num_params;
    %
    %> Number of veiling variables \mathbf{v} of the system.
    %
    m_num_veils;
    %
  end
  %
  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Class constructor for a system:
    %>
    %> \param t_data          The system data.
    %> \param t_num_free_dofs The number of free dofs of the system.
    %> \param t_num_spec_dofs The number of specified dofs of the system.
    %> \param t_num_params    The number of parameters \mathbf{x} of the system.
    %> \param t_num_veils     The number of veiling variables \mathbf{v} of the system.
    %>
    %> \return The system.
    %
    function this = System( t_data, t_num_free_dofs, t_num_spec_dofs, t_num_params, t_num_veils )
      this.m_data          = t_data;
      this.m_num_free_dofs = t_num_free_dofs;
      this.m_num_spec_dofs = t_num_spec_dofs;
      this.m_num_params    = t_num_params;
      this.m_num_veils     = t_num_veils;
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
    %> Get the system data field.
    %>
    %> \param field The system data field.
    %>
    %> \return The system data field.
    %
    function out = get_data_field( this, field )
      out = this.m_data.(field);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Set the system data field.
    %>
    %> \param field The system data field.
    %> \param value The system data field value.
    %
    function set_data_field( this, field, value )
      this.m_data.(field) = value;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Get the number of free dofs of the system.
    %>
    %> \return The number of free dofs of the system.
    %
    function out = get_num_free_dofs( this )
      out = this.m_num_free_dofs;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Set the number of free dofs of the system.
    %>
    %> \param t_num_free_dofs The number of free dofs of the system.
    %
    function set_num_free_dofs( this, t_num_free_dofs )
      this.m_num_free_dofs = t_num_free_dofs;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Get the number of specified dofs of the system.
    %>
    %> \return The number of specified dofs of the system.
    %
    function out = get_num_spec_dofs( this )
      out = this.m_num_spec_dofs;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Set the number of specified dofs of the system.
    %>
    %> \param t_num_spec_dofs The number of specified dofs of the system.
    %
    function set_num_spec_dofs( this, t_num_spec_dofs )
      this.m_num_spec_dofs = t_num_spec_dofs;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Get the number of parameters \mathbf{x} of the system.
    %>
    %> \return The number of parameters \mathbf{x} of the system.
    %
    function out = get_num_params( this )
      out = this.m_num_params;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Set the number of parameters \mathbf{x} of the system.
    %>
    %> \param t_num_params The number of parameters \mathbf{x} of the system.
    %
    function set_num_params( this, t_num_params )
      this.m_num_params = t_num_params;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Get the number of veiling variables \mathbf{c} of the system.
    %>
    %> \return The number of veiling variables \mathbf{c} of the system.
    %
    function out = get_num_veils( this )
      out = this.m_num_veils;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Set the number of veiling variables \mathbf{c} of the system.
    %>
    %> \param t_num_params The number of veiling variables \mathbf{c} of the system.
    %
    function set_num_veils( this, t_num_params )
      this.m_num_veils = t_num_veils;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the system deformation vector \f$ \mathbf{d} \f$ as:
    %>
    %> \f[
    %>   \mathbf{d} = \left[
    %>     \mathbf{d}_{f}, \mathbf{d}_{s}
    %>   \right]^{\top}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The system deformation vector \f$ \mathbf{d} \f$.
    %
    function out = compute_d( this, x, v, varargin )

      % Compute specified deformations
      d_s = this.d_s(x, v);

      % Compute free deformations
      if nargin == 3
        d_f = this.K_ff(x, v)\(this.f_f(x, v)-this.K_fs(x, v)*d_s); ...
      elseif nargin == 5
        [d_f, ~] = lsqr(this.K_ff(x, v), this.f_f(x, v)-this.K_fs(x, v)*d_s, varargin{1}, varargin{2});
      else
        error('TrussMe.System.compute_d(...): Wrong number of arguments.');
      end

      % Permute the result
      out = [d_f; d_s];
      out = out(this.unperm());
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
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The system deformation vector \f$ \mathbf{d}_{f} \f$.
    %
    function out = compute_d_f( this, x, v, varargin )
      if nargin == 3
        out = this.K_ff(x, v)\(this.f_f(x, v)-this.K_fs(x, v)*this.d_s(x, v)); ...
      elseif nargin == 5
        [out, ~] = lsqr(this.K_ff(x, v), this.f_f(x, v)-this.K_fs(x, v)*this.d_s(x, v), varargin{1}, varargin{2});
      else
        error('TrussMe.System.compute_d_f(...): Wrong number of arguments.');
      end
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the deformation vector \f$ \mathbf{d} \f$ with
    %> respect to \f$ \mathbf{x} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jd}_{\mathbf{x}} = \left[
    %>     \mathbf{Jd}_{f\mathbf{x}}, \mathbf{Jd}_{s\mathbf{x}}
    %>   \right]^{\top}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the deformation vector \f$ \mathbf{Jd}_{\mathbf{x}} \f$.
    %
    function out = compute_Jd_x( this, x, v, varargin )
      out = [ ...
        this.Jd_f_x(x, v, varargin{:}); ...
        this.compute_Jd_s_x(x, v) ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{d}_{f} \f$ with respect to \f$ \mathbf{x} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jd}_{f\mathbf{x}} =
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{Jd}_{f\mathbf{x}} \f$.
    %
    function out = compute_Jd_f_x( this, x, v, varargin )
      out = [];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the deformation vector \f$ \mathbf{d} \f$ with
    %> respect to \f$ \mathbf{v} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jd}_{\mathbf{v}} = \left[
    %>     \mathbf{Jd}_{f\mathbf{v}}, \mathbf{Jd}_{s\mathbf{v}}
    %>   \right]^{\top}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the deformation vector \f$ \mathbf{Jd}_{\mathbf{v}} \f$.
    %
    function out = compute_Jd_v( this, x, v, varargin )
      out = [ ...
        this.Jd_f_v(x, v, varargin{:}); ...
        this.compute_Jd_s_v(x, v) ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{d}_{f} \f$ with respect to \f$ \mathbf{v} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jd}_{f\mathbf{v}} =
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{Jd}_{f\mathbf{v}} \f$.
    %
    function out = compute_Jd_f_v( this, x, v, varargin )
      out = [];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the system force vector \f$ \mathbf{f} \f$ as:
    %>
    %> \f[
    %>   \mathbf{f} = \left[
    %>     \mathbf{f}_{f}, \mathbf{f}_{s}
    %>   \right]^{\top}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The system force vector \f$ \mathbf{f} \f$.
    %
    function out = compute_f( this, x, v, varargin )
      out = [ ...
        this.f_f(x, v); ...
        this.K_sf(x, v)*this.compute_d_f(x, v, varargin{:}) + ...
        this.K_ss(x, v)*this.d_s(x, v) - this.f_r(x, v); ...
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
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The system force vector \f$ \mathbf{f}_{s} \f$.
    %
    function out = compute_f_s( this, x, v, varargin )
      out = this.K_sf(x, v)*this.compute_d_f(x, v, varargin{:}) + ...
            this.K_ss(x, v)*this.d_s(x, v) - this.f_r(x, v);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the system force vector \f$ \mathbf{f} \f$ with
    %> respect to \f$ \mathbf{x} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jf}_{\mathbf{x}} = \left[
    %>     \mathbf{Jf}_{f\mathbf{x}}, \mathbf{Jf}_{s\mathbf{x}}
    %>   \right]^{\top}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the system force vector
    %> \f$ \mathbf{Jf}_{\mathbf{x}} \f$.
    %
    function out = compute_Jf_x( this, x, v, varargin )
      out = [ ...
        this.Jf_f_x(x, v); ...
        this.compute_Jf_s_x(x, v, varargin{:}) ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the system force vector of specified dofs
    %> \f$ \mathbf{f}_{s} \f$ with respect to \f$ \mathbf{x} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jf}_{s\mathbf{x}} = \mathbf{TK}_{sf\mathbf{x}} \mathbf{d}_{f} +
    %>     \mathbf{K}_{sf} \mathbf{Jd}_{f\mathbf{x}} +
    %>     \mathbf{TK}_{ss\mathbf{x}} \mathbf{d}_{s} +
    %>     \mathbf{K}_{ss} \mathbf{Jd}_{s\mathbf{x}} - \mathbf{Jf}_{r\mathbf{x}}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{Jf}_{s\mathbf{x}} \f$.
    %
    function out = compute_Jf_s_x( this, x, v, varargin )
      out = zeros(this.m_num_free_dofs, this.m_num_params);
      out = this.TK_sf_x(x, v)*this.compute_d_f(x, v, varargin{:}) + ... % FIXME
            this.K_sf(x, v)*this.compute_Jd_f_x(x, v, varargin{:}) + ...
            this.TK_ss_x(x, v)*this.d_s(x, v) + ...
            this.K_ss(x, v)*this.Jd_s_x(x, v) - this.Jf_r_x(x, v);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the system force vector \f$ \mathbf{f} \f$ with
    %> respect to \f$ \mathbf{v} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jf}_{\mathbf{v}} = \left[
    %>     \mathbf{Jf}_{f\mathbf{v}}, \mathbf{Jf}_{s\mathbf{v}}
    %>   \right]^{\top}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the system force vector
    %> \f$ \mathbf{Jf}_{\mathbf{v}} \f$.
    %
    function out = compute_Jf_v( this, x, v, varargin )
      out = [ ...
        this.Jf_f_v(x, v); ...
        this.compute_Jf_s_v(x, v, varargin{:}) ...
      ];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Compute the Jacobian of the system force vector of specified dofs
    %> \f$ \mathbf{f}_{s} \f$ with respect to \f$ \mathbf{v} \f$ as:
    %>
    %> \f[
    %>   \mathbf{Jf}_{s\mathbf{v}} = \mathbf{TK}_{sf\mathbf{v}} \mathbf{d}_{f} +
    %>     \mathbf{K}_{sf} \mathbf{Jd}_{f\mathbf{v}} +
    %>     \mathbf{TK}_{ss\mathbf{v}} \mathbf{d}_{s} +
    %>     \mathbf{K}_{ss} \mathbf{Jd}_{s\mathbf{v}} - \mathbf{Jf}_{r\mathbf{v}}
    %> \f]
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return The Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{Jf}_{s\mathbf{v}} \f$.
    %
    function out = compute_Jf_s_v( this, x, v, varargin )
      out = this.TK_sf_v(x, v)*this.compute_d_f(x, v, varargin{:}) + ... % FIXME
            this.K_sf(x, v)*this.compute_Jd_f_v(x, v, varargin{:}) + ...
            this.TK_ss_v(x, v)*this.d_s(x, v) + ...
            this.K_ss(x, v)*this.Jd_s_v(x, v) - this.Jf_r_v(x, v);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Internal structure sanity check.
    %>
    %> \param x   States \f$ \mathbf{x} \f$.
    %> \param v   Veils \f$ \mathbf{v} \f$.
    %> \param tol [optional] Tolerance for the iterative solver.
    %> \param itr [optional] Maximum number of solver iterations.
    %>
    %> \return An error is thrown if the sizes are not correct.
    %
    function sanity_check( this, x, v, varargin )

      CMD = 'TrussMe.System.sanity_check(...): ';

      % Evaluate stiffness matrices
      K_ff = this.K_ff(x, v);
      K_sf = this.K_sf(x, v);
      K_fs = this.K_fs(x, v);
      K_ss = this.K_ss(x, v);
      K    = this.K(x, v);

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
      d_fc = this.compute_d_f(x, v, varargin{:});
      d_c = this.compute_d(x, v, varargin{:});

      % Compute force vectors
      f_sc = this.compute_f_s(x, v);
      f_c  = this.compute_f(x, v);

      % Check sizes
      assert(size(K_ff, 1) == size(K_ff, 2), [CMD, 'K_ff is not square.']);
      assert(size(K_ss, 1) == size(K_ss, 2), [CMD, 'K_ss is not square.']);
      assert(size(K, 1)    == size(K, 2),    [CMD, 'K is not square.']);

      assert(size(K_ff, 1) + size(K_sf, 1) == size(K, 1) && ...
             size(K_ff, 2) + size(K_fs, 2) == size(K, 2), ...
        [CMD, 'K_ff, K_sf and K are not compatible.']);
      assert(size(K_ss, 1) + size(K_fs, 1) == size(K, 1) && ...
             size(K_ss, 2) + size(K_sf, 2) == size(K, 2), ...
        [CMD, 'K_ss, K_fs and K are not compatible.']);

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
    %> Evaluate the tensor of the system stiffness matrix \f$ \mathbf{K} \f$
    %> with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{\mathbf{x}} \f$.
    %
    TK_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the free-free dofs stiffness matrix
    %> \f$ \mathbf{K}_{ff} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{ff\mathbf{x}} \f$.
    %
    TK_ff_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the free-specified dofs stiffness matrix
    %> \f$ \mathbf{K}_{fs} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{fs\mathbf{x}} \f$.
    %
    TK_fs_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the specified-free dofs stiffness matrix
    %> \f$ \mathbf{K}_{sf} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{sf\mathbf{x}} \f$.
    %
    TK_sf_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the specified-specified dofs stiffness matrix
    %> \f$ \mathbf{K}_{ss} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{ss\mathbf{x}} \f$.
    %
    TK_ss_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the system stiffness matrix \f$ \mathbf{K} \f$
    %> with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{\mathbf{v}} \f$.
    %
    TK_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the free-free dofs stiffness matrix
    %> \f$ \mathbf{K}_{ff} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{ff\mathbf{v}} \f$.
    %
    TK_ff_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the free-specified dofs stiffness matrix
    %> \f$ \mathbf{K}_{fs} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{fs\mathbf{v}} \f$.
    %
    TK_fs_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the specified-free dofs stiffness matrix
    %> \f$ \mathbf{K}_{sf} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{sf\mathbf{v}} \f$.
    %
    TK_sf_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the tensor of the specified-specified dofs stiffness matrix
    %> \f$ \mathbf{K}_{ss} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The tensor \f$ \mathbf{TK}_{ss\mathbf{v}} \f$.
    %
    TK_ss_v( this, x, v )
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
    %> Evaluate the Jacobian of the deformation vector \f$ \mathbf{d} \f$ with
    %> respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{Jd}_{\mathbf{x}} \f$.
    %
    Jd_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{d}_{f} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{Jd}_{f\mathbf{x}} \f$.
    %
    Jd_f_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the specified dofs deformation vector
    %> \f$ \mathbf{d}_{s} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{Jd}_{s\mathbf{x}} \f$.
    %
    Jd_s_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the deformation vector \f$ \mathbf{d} \f$ with
    %> respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{Jd}_{\mathbf{v}} \f$.
    %
    Jd_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the free dofs deformation vector
    %> \f$ \mathbf{d}_{f} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{Jd}_{f\mathbf{v}} \f$.
    %
    Jd_f_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the specified dofs deformation vector
    %> \f$ \mathbf{d}_{s} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system deformation vector \f$ \mathbf{Jd}_{s\mathbf{v}} \f$.
    %
    Jd_s_v( this, x, v )
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
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system force vector \f$ \mathbf{f} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{Jf}_{\mathbf{x}} \f$.
    %
    Jf_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system force vector of free dofs \f$ \mathbf{f}_{f} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{Jf}_{f\mathbf{x}} \f$.
    %
    Jf_f_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system force vector of specified dofs \f$ \mathbf{f}_{s} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{Jf}_{s\mathbf{x}} \f$.
    %
    Jf_s_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system remainder force vector
    %> \f$ \mathbf{f}_{r} \f$ with respect to \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system remainder force vector \f$ \mathbf{Jf}_{r\mathbf{x}} \f$.
    %
    Jf_r_x( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system force vector \f$ \mathbf{f}
    %> \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{Jf}_{\mathbf{v}} \f$.
    %
    Jf_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system force vector of free dofs
    %> \f$ \mathbf{f}_{f} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{Jf}_{f\mathbf{v}} \f$.
    %
    Jf_f_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system force vector of specified dofs
    %> \f$ \mathbf{f}_{s} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system force vector \f$ \mathbf{Jf}_{s\mathbf{v}} \f$.
    %
    Jf_s_v( this, x, v )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the system remainder force vector
    %> \f$ \mathbf{f}_{r} \f$ with respect to \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %> \param v Veils \f$ \mathbf{v} \f$.
    %>
    %> \return The system remainder force vector \f$ \mathbf{Jf}_{r\mathbf{v}} \f$.
    %
    Jf_r_v( this, x, v )
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
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the veils \f$ \mathbf{v} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %>
    %> \return The Veils \f$ \mathbf{v} \f$.
    %
    v( this, x )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %> Evaluate the Jacobian of the veils \f$ \mathbf{v} \f$ with respect to
    %> \f$ \mathbf{x} \f$.
    %>
    %> \param x States \f$ \mathbf{x} \f$.
    %>
    %> \return The Veils \f$ \mathbf{v} \f$.
    %
    Jv_x( this, x )
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
  %
end

% That's All Folks!
