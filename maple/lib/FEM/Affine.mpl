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

export IsFRAME := proc(
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

export GenerateFrame := proc(
  p_1::{list, Vector, VECTOR, POINT},
  p_2::{list, Vector, VECTOR, POINT},
  vec::{list, Vector, VECTOR},
  $)::FRAME;

  description "Generate a reference frame matrix from three points or vectors "
    "<p_1>, <p_2> and vector <vec> ortogonal to XY-plane <vec>";

  local p_1_tmp, p_2_tmp, vec_tmp, e_x, e_y, e_z;

  if type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 3) then
    p_1_tmp := p_1;
  elif type(p_1, list) and evalb(nops(p_1) = 3) then
    p_1_tmp := convert(p_1, Vector);
  elif type(p_1, VECTOR) or type(p_1, POINT) then
    p_1_tmp := p_1[1..3];
  else
    error("<p_1> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  if type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 3) then
    p_2_tmp := p_2;
  elif type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp := convert(p_2, Vector);
  elif type(p_2, VECTOR) or type(p_2, POINT) then
    p_2_tmp := p_2[1..3];
  else
    error("<p_2> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  if type(vec, Vector) and evalb(LinearAlgebra:-Dimension(vec) = 3) then
    vec_tmp := vec;
  elif type(vec, list) and evalb(nops(vec) = 3) then
    vec_tmp := convert(vec, Vector);
  elif type(vec, VECTOR) or type(vec, POINT) then
    vec_tmp := vec[1..3];
  else
    error("<vec> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  e_x := p_2_tmp - p_1_tmp;
  e_x := e_x / TrussMe:-FEM:-Norm2(e_x);

  e_y := LinearAlgebra:-CrossProduct(vec_tmp, e_x);
  e_y := e_y / TrussMe:-FEM:-Norm2(e_y);

  e_z := LinearAlgebra:-CrossProduct(e_x, e_y);
  e_z := e_z / TrussMe:-FEM:-Norm2(e_z);

  return <<e_x,     0>|
          <e_y,     0>|
          <e_z,     0>|
          <p_1_tmp, 1>>;
end proc: # GenerateFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateGenericFrame := proc(
  label::string := "",
  {
  e := "e",
  p := "p",
  x := "x",
  y := "y",
  z := "z",
  s := "__"
  }, $)::FRAME;

  description "Generate a generic reference frame matrix from a string label "
    "<label>. Optional arguments <e>, <p>, <x>, <y>, <z> and <s> are used to "
    "customize the output, e.g., 'exx := e||s||x||x'.";

  local exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz, pxx, pyy, pzz;

  exx := e||s||x||x; eyx := e||s||y||x; ezx := e||s||z||x;
  exy := e||s||x||y; eyy := e||s||y||y; ezy := e||s||z||y;
  exz := e||s||x||z; eyz := e||s||y||z; ezz := e||s||z||z;
  pxx := p||s||x||x; pyy := p||s||y||y; pzz := p||s||z||z;

  if evalb(label <> "") then
    exx := cat(exx, s||label); eyx := cat(eyx, s||label); ezx := cat(ezx, s||label);
    exy := cat(exy, s||label); eyy := cat(eyy, s||label); ezy := cat(ezy, s||label);
    exz := cat(exz, s||label); eyz := cat(eyz, s||label); ezz := cat(ezz, s||label);
    pxx := cat(pxx, s||label); pyy := cat(pyy, s||label); pzz := cat(pzz, s||label);
  end if;

  return <<convert(exx, symbol), convert(exy, symbol), convert(exz, symbol), 0>|
          <convert(eyx, symbol), convert(eyy, symbol), convert(eyz, symbol), 0>|
          <convert(ezx, symbol), convert(ezy, symbol), convert(ezz, symbol), 0>|
          <convert(pxx, symbol), convert(pyy, symbol), convert(pzz, symbol), 1>>;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export InverseFrame := proc(
  RF::FRAME,
  $)::FRAME;

  description "Inverse affine transformation matrix <RF>.";

  LinearAlgebra:-Transpose(RF[1..3, 1..3]);
  return <<% | -%.RF[1..3, 4]>,
          <0 | 0 | 0 | 1>>;
end proc: # InverseFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Translate := proc(
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

export  Translation := proc(
  RF::FRAME,
  $)::Vector;

  description "Extract the translation vector of the reference frame <RF>.";

  return RF[1..3, 4];
end proc: # Translation

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Rotate := proc(
  axis::{symbol, string},
  angle::algebraic,
  $)::FRAME;

  description "Affine transformation matrix corresponding to the rotation "
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
    error("<axis> must be a string or a symbol.");
  end if;
end proc: # Rotate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Rotation := proc(
  RF::FRAME,
  $)::Matrix;

  description "Extract the rotation matrix of the reference frame <RF>.";

  return RF[1..3, 1..3];
end proc: # Rotation

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsVECTOR := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of VECTOR type.";

  return type(var, Vector) and
    evalb(LinearAlgebra:-Dimension(var) = 4) and evalb(var[4] = 0);
end proc: # IsVECTOR

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsPOINT := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of POINT type.";

  return type(var, Vector) and
    evalb(LinearAlgebra:-Dimension(var) = 4) and evalb(var[4] = 1);
end proc: # IsPOINT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Origin := proc(
  RF::FRAME,
  $)::POINT;

  description "Extract the origin point of the reference frame <RF>.";

  return convert(<RF[1..3, 4], 1>, Vector);
end proc: # Origin

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CompX := proc(
  x::{VECTOR, POINT},
  $)::algebraic;

  description "Extract the x-axis component of the vector or point <x>.";

  return x[1];
end proc: # CompX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CompY := proc(
  x::{VECTOR, POINT},
  $)::algebraic;

  description "Extract the y-axis component of the vector or point <x>.";

  return x[2];
end proc: # CompY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CompZ := proc(
  x::{VECTOR, POINT},
  $)::algebraic;

  description "Extract the z-axis component of the vector or point <x>.";

  return x[3];
end proc: # CompZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CompXYZ := proc(
  x::{VECTOR, POINT},
  $)::algebraic, algebraic, algebraic;

  description "Extract the x, y and z-axis components of the vector or point "
    "<x>.";

  return x[1], x[2], x[3];
end proc: # CompXYZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecX := proc(
  RF::FRAME := Matrix(4, shape = identity),
  $)::VECTOR;

  description "Extract the x-axis unit vector of the reference frame <RF>.";

  return convert(<RF[1..3, 1], 0>, Vector);
end proc: # UvecX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecY := proc(
  RF::FRAME,
  $)::VECTOR;

  description "Extract the y-axis unit vector of the reference frame <RF>.";

  return convert(<RF[1..3, 2], 0>, Vector);
end proc: # UvecY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecZ := proc(
  RF::FRAME := Matrix(4, shape = identity),
  $)::VECTOR;

  description "Extract the z-axis unit vector of the reference frame <RF>.";

  return convert(<RF[1..3, 3], 0>, Vector);
end proc: # UvecZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecXYZ := proc(
  RF::FRAME := Matrix(4, shape = identity),
  $)::VECTOR, VECTOR, VECTOR;

  description "Extract the x, y and z-axis unit vectors of the reference "
    "frame <RF>.";

  return convert(<RF[1..3, 1], 0>, Vector),
         convert(<RF[1..3, 2], 0>, Vector),
         convert(<RF[1..3, 3], 0>, Vector);
end proc: # UvecXYZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Project := proc(
  x::{VECTOR, POINT},
  RF_ini::FRAME,
  RF_end::FRAME,
  $)::{VECTOR, POINT};

  description "Project the vector or point <x>, from reference frame <RF_ini> "
    "to reference frame <RF_end>.";

  # Try to compare reference frames
  try
    # FIXME: problems with floats (floats not handled error)
    evalb~(evala(TrussMe:-FEM:-Simplify(RF_end) =~ TrussMe:-FEM:-Simplify(RF_ini)));
  catch:
    evalb~(RF_end =~ RF_ini);
  end try;

  # Return the projection
  if has(%, false) then
    return TrussMe:-FEM:-Simplify(TrussMe:-FEM:-InverseFrame(RF_end).RF_ini.x);
  else
    return x;
  end if;
end proc: # Project

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
