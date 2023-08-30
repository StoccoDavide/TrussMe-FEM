# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                    ____          _                                          #
#                   / ___|___   __| | ___  __ _  ___ _ __                     #
#                  | |   / _ \ / _` |/ _ \/ _` |/ _ \ '_ \                    #
#                  | |__| (_) | (_| |  __/ (_| |  __/ | | |                   #
#                   \____\___/ \__,_|\___|\__, |\___|_| |_|                   #
#                                         |___/                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DoDiff := proc(
  # _passed
  )::anything;

  description "Differentiate an expression with respect to a function.";

  local i, tmp1, tmp2, out;

  tmp2 := diff(subs(f = tmp1, out), tmp1);
  subs(tmp1 = f, convert(tmp2, D));
  out := _passed[1];
  for i from 2 to _npassed do
    tmp2 := diff(subs(_passed[i] = tmp1, out), tmp1);
    out  := subs(tmp1 = _passed[i], convert(tmp2, D));
  end do;
  return out;
end proc: # DoDiff

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DoGradient := proc(
  fnc::algebraic,
  lst::{Vector, list},
  $)::Vector;

  description "Differentiate a scalar expression <fnc> with respect to a "
    "list of functions <lst>.";

  local i, n, out;

  # Extract dimensions
  if type(lst, Vector) then
    n := LinearAlgebra:-Dimension(lst);
  else
    n := nops(lst);
  end if;

  # Differentiate
  out := Vector(n);
  for i from 1 to n do
    out[i] := TrussMe_FEM:-DoDiff(fnc, lst[i]):
  end do;
  return out;
end proc: # DoGradient

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DoHessian := proc(
  fnc::algebraic,
  lst::{Vector, list},
  $)::Matrix;

  description "Differentiate a vector of expressions (gradient) <fnc> with "
    "respect to a list <lst> of functions.";

  local i, j, n, out;

  # Extract dimensions
  if type(lst, Vector) then
    n := LinearAlgebra:-Dimension(lst);
  else
    n := nops(lst);
  end if;

  # Differentiate
  out := Matrix(n, n);
  for i from 1 to n do
    for j from i to n do
      out[i, j] := TrussMe_FEM:-DoDiff(fnc, lst[i], lst[j]);
      out[j, i] := out[i, j];
    end do;
  end do;
  return out;
end proc: # DoHessian

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DoJacobian := proc(
  fnc::Vector,
  lst::{Vector, list},
  $)::Matrix;

  description "Differentiate a vector of expressions <fnc> with respect to a "
    "list <lst> of functions.";

  local i, j, m, n, out;

  # Extract dimensions
  m := LinearAlgebra:-Dimension(fnc);
  if type(lst, Vector) then
    n := LinearAlgebra:-Dimension(lst);
  else
    n := nops(lst);
  end if;

  # Differentiate
  out := Matrix(m, n);
  for i from 1 to m do
    for j from 1 to n do
      out[i, j] := TrussMe_FEM:-DoDiff(fnc[i], lst[j]);
    end do;
  end do;
  return out;
end proc: # DoJacobian

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DoTensor := proc(
  mat::Matrix,
  lst::{Vector, list},
  $)::Array;

  description "Differentiate a matrix <mat> with respect to a list of "
    "functions <lst>.";

  local i, j, k, m, n, l, out;

  # Extract dimensions
  m := LinearAlgebra:-RowDimension(mat);
  n := LinearAlgebra:-ColumnDimension(mat);
  if type(lst, Vector) then
    l := LinearAlgebra:-Dimension(lst);
  else
    l := nops(lst);
  end if;

  # Differentiate
  out := Array(1..m, 1..n, 1..l);
  for i from 1 to m do
    for j from 1 to n do
      for k from 1 to l do
        out[i, j, k] := TrussMe_FEM:-DoDiff(mat[i, j], lst[k]);
      end do;
    end do;
  end do;
  return out;
end proc: # DoTensor

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
