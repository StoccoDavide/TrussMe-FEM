# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           _     __  __ _                                    #
#                          / \   / _|/ _(_)_ __   ___                         #
#                         / _ \ | |_| |_| | '_ \ / _ \                        #
#                        / ___ \|  _|  _| | | | |  __/                        #
#                       /_/   \_\_| |_| |_|_| |_|\___|                        #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsFRAME::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of FRAME type.";

  return type(var, Matrix) and
    evalb(LinearAlgebra:-RowDimension(var) = 4) and
    evalb(LinearAlgebra:-ColumnDimension(var) = 4) and
    evalb(var[4, 1] = 0) and evalb(var[4, 2] = 0) and
    evalb(var[4, 3] = 0) and evalb(var[4, 4] = 1);
end proc: # IsFRAME

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export InverseFrame::static := proc(
  RF::FRAME,
  $)::FRAME;

  description "Inverse of the affine transformation matrix <RF>.";

  LinearAlgebra:-Transpose(RF[1..3, 1..3]);
  return <<% | -%.RF[1..3, 4]>,
          <0 | 0 | 0 | 1>>;
end proc: # InverseFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Translate::static := proc(
  x::algebraic,
  y::algebraic,
  z::algebraic,
  $)::FRAME;

  description "Affine transformation matrix corresponding to a translation "
    "<x, y, z>.";

  return <<1, 0, 0, 0>|
          <0, 1, 0, 0>|
          <0, 0, 1, 0>|
          <x, y, z, 1>>;
end proc: # Translate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Rotate::static := proc(
  axis::{symbol, string},
  angle::algebraic,
  $)::FRAME;

  description "Affinet ransformation matrix corresponding to the rotation "
    "<angle> around the given <axis>.";

  if evalb(axis = 'X') or evalb(axis = "X") then
    return <<1, 0,           0,          0>|
            <0, cos(angle),  sin(angle), 0>|
            <0, -sin(angle), cos(angle), 0>|
            <0, 0,           0,          1>>;
  elif evalb(axis = 'Y') or evalb(axis = "Y") then
    return <<cos(angle), 0, -sin(angle), 0>|
            <0,          1, 0,           0>|
            <sin(angle), 0, cos(angle),  0>|
            <0,          0, 0,           1>>;
  elif evalb(axis = 'Z') or evalb(axis = "Z") then
    return <<cos(angle),  sin(angle), 0, 0>|
            <-sin(angle), cos(angle), 0, 0>|
            <0,           0,          1, 0>|
            <0,           0,          0, 1>>;
  else
    error "invalid axis detected.";
  end if;
end proc: # Rotate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsVECTOR::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of VECTOR type.";

  return type(var, Vector) and
    evalb(LinearAlgebra:-Dimension(var) = 4) and evalb(var[4] = 0);
end proc: # IsVECTOR

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsPOINT::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of POINT type.";

  return type(var, Vector) and
    evalb(LinearAlgebra:-Dimension(var) = 4) and evalb(var[4] = 1);
end proc: # IsPOINT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Origin::static := proc(
  RF::FRAME,
  $)::VECTOR;

  description "Extract the origin point of the reference frame <RF>.";

  return <RF[1..3, 4], 1>;
end proc: # Origin

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Uvec::static := proc(
  axis::{symbol, string},
  RF::FRAME := Matrix(4, shape = identity),
  $)::VECTOR;

  description "Extract the unit vector of the reference frame <RF> along the "
    "given <axis>.";

  if evalb(axis = 'X') or evalb(axis = "X") then
    return TrussMe_FEM:-UvecX(RF);
  elif evalb(axis = 'Y') or evalb(axis = "Y") then
    return TrussMe_FEM:-UvecY(RF);
  elif evalb(axis = 'Z') or evalb(axis = "Z") then
    return TrussMe_FEM:-UvecZ(RF);
  else
    error "invalid axis detected.";
  end if;
end proc: # Uvec

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecX::static := proc(
  RF::FRAME := Matrix(4, shape = identity),
  $)::VECTOR;

  description "Extract the x-axis unit vector of the reference frame <RF>.";

  return <RF[1..3, 1], 0>;
end proc: # UvecX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecY::static := proc(
  RF::FRAME,
  $)::VECTOR;

  description "Extract the y-axis unit vector of the reference frame <RF>.";

  return <RF[1..3, 2], 0>;
end proc: # UvecY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecZ::static := proc(
  RF::FRAME := Matrix(4, shape = identity),
  $)::VECTOR;

  description "Extract the z-axis unit vector of the reference frame <RF>.";

  return <RF[1..3, 3], 0>;
end proc: # UvecZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Project::static := proc(
  x::{VECTOR, POINT},
  RF_ini::FRAME,
  RF_end::FRAME,
  $)::{VECTOR, POINT};

  description "Project the vector or point <x>, from reference frame <RF_ini> "
    "to reference frame <RF_end>.";

  # Try to compare reference frames
  try
    # FIXME: problems with floats (floats not handled error)
    evalb~(evala(TrussMe_FEM:-Simplify(RF_end) =~ TrussMe_FEM:-Simplify(RF_ini)));
  catch:
    evalb~(RF_end =~ RF_ini);
  end try;

  # Return the projection
  if has(%, false) then
    return TrussMe_FEM:-Simplify(TrussMe_FEM:-InverseFrame(RF_end).RF_ini.x);
  else
    return x;
  end if;
end proc: # Project

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
