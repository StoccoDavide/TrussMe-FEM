# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                ____  _                   _                                  #
#               / ___|| |_ _ __ _   _  ___| |_ _   _ _ __ ___                 #
#               \___ \| __| '__| | | |/ __| __| | | | '__/ _ \                #
#                ___) | |_| |  | |_| | (__| |_| |_| | | |  __/                #
#               |____/ \__|_|   \__,_|\___|\__|\__,_|_|  \___|                #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsNODE::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of NODE type.";

  return type(var, table) and evalb(var["type"] = NODE);
end proc: # IsNODE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSUPPORT::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of SUPPORT type.";

  return type(var, NODE) and has(map(has, var["dofs"], 0), true);
end proc: # IsSUPPORT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsDOFS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of DOFS type.";

  return type(var, list) and evalb(nops(var) = 6) and
    not has(map(has, var, {0, 1}), false);
end proc: # IsDOFS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeNode::static := proc(
  name::string,
  coords::{list, POINT, VECTOR},
  RF::FRAME := Matrix(4, shape = identity),
  {
  dofs::DOFS := [1, 1, 1, 1, 1, 1]
  }, $)::NODE;

  description "Create a node with name <name> at coordinates <coords> in the "
    "reference frame <RF>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained to the ground "
    "in the direction given from <RF>.";

  local coords_tmp;

  if type(coords, list) and evalb(nops(coords) = 3) then
    coords_tmp := coords;
  elif type(coords, POINT) or type(coords, VECTOR) then
    coords_tmp := [TrussMe_FEM:-CompXYZ(coords)];
  else
    error("<coords> must be a list of 3 elements, a POINT or a VECTOR.");
  end if;

  return table([
    "type"                = NODE,
    "name"                = name,
    "id"                  = TrussMe_FEM:-GenerateId(),
    "frame"               = RF,
    "coordinates"         = coords_tmp,
    "dofs"                = dofs,
    "displacements"       = [],
    "output_reactions"    = [], # Output loads in the node frame
    "output_deformations" = []  # Output displacements in the node frame
  ]);
end proc: # MakeNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCompliantNode := proc(
  name::string,
  coords::{list, POINT, VECTOR},
  RF::FRAME := Matrix(4, shape = identity),
  {
  dofs::DOFS := [1, 1, 1, 1, 1, 1],
  K::{algebraic, list(algebraic)} := 0,
  T::{algebraic, list(algebraic)} := 0
  }, $)::NODE, ELEMENT, NODE;

  description "Create a node with name <name> at coordinates <coords> in the "
    "reference frame <RF>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained to the ground "
    "in the direction given from <RF>. The node is also connected to a "
    "compliant spring element with traslational stiffness <K> and torsional "
    "stiffness <T>.";

  local fixed, spring, compliant;

  fixed     := TrussMe_FEM:-MakeNode(name, coords, RF, parse("dofs") = dofs);
  compliant := TrussMe_FEM:-MakeNode(cat(name, "_compliant"), coords, RF);
  spring    := TrussMe_FEM:-MakeSpring(
    cat(name, "_spring"), RF, fixed, compliant, parse("K") = K, parse("T") = T
  );
  return fixed, spring, compliant;
end proc: # MakeCompliantNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSTIFFNESS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of STIFFNESS type.";

  return type(var, Matrix) and
    evalb(LinearAlgebra:-RowDimension(var) = 12) and
    evalb(LinearAlgebra:-ColumnDimension(var) = 12);
end proc: # IsSTIFFNESS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetSpringStiffness::static := proc(
  K::{algebraic, list(algebraic)},
  T::{algebraic, list(algebraic)},
  $)::STIFFNESS;

  description "Get the stiffness matrix of a spring given the translational "
    "stiffnesses <K> or <K_x, K_y, K_z>, and the torsional stiffnesses <T> or "
    "<T_x, T_y, T_z>.";

  local K_x, K_y, K_z, T_x, T_y, T_z;

  if type(K, list(algebraic)) and evalb(nops(K) = 3) then
    K_x := K[1]; K_y := K[2]; K_z := K[3];
  elif type(K, algebraic) then
    K_x := K; K_y := K; K_z := K;
  else
    error "<K> must be a list of 3 elements or a single element.";
  end if;

  if type(T, list(algebraic)) and evalb(nops(T) = 3) then
    T_x := T[1]; T_y := T[2]; T_z := T[3];
  elif type(T, algebraic) then
    T_x := T; T_y := T; T_z := T;
  else
    error "<T> must be a list of 3 elements or a single element.";
  end if;

  return Matrix(
    <<K_x, 0, 0, 0, 0, 0, -K_x, 0, 0, 0, 0, 0>|
     <0, K_y, 0, 0, 0, 0, 0, -K_y, 0, 0, 0, 0>|
     <0, 0, K_z, 0, 0, 0, 0, 0, -K_z, 0, 0, 0>|
     <0, 0, 0, T_x, 0, 0, 0, 0, 0, -T_x, 0, 0>|
     <0, 0, 0, 0, T_y, 0, 0, 0, 0, 0, -T_y, 0>|
     <0, 0, 0, 0, 0, T_z, 0, 0, 0, 0, 0, -T_z>|
     <-K_x, 0, 0, 0, 0, 0, K_x, 0, 0, 0, 0, 0>|
     <0, -K_y, 0, 0, 0, 0, 0, K_y, 0, 0, 0, 0>|
     <0, 0, -K_z, 0, 0, 0, 0, 0, K_z, 0, 0, 0>|
     <0, 0, 0, -T_x, 0, 0, 0, 0, 0, T_x, 0, 0>|
     <0, 0, 0, 0, -T_y, 0, 0, 0, 0, 0, T_y, 0>|
     <0, 0, 0, 0, 0, -T_z, 0, 0, 0, 0, 0, T_z>>,
    storage = sparse);
end proc: # GetSpringStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetRodStiffness::static := proc(
  A::algebraic,
  E::algebraic,
  G::algebraic,
  L::algebraic,
  I_x::algebraic,
  I_y::algebraic,
  I_z::algebraic,
  $)::STIFFNESS;

  description "Get the stiffness matrix of a rod (only axial-stiffness) "
    "given the cross-section area <A>, the Young's modulus <E>, the shear "
    "modulus <G>, the length <L>, and the cross-section inertia "
    "<I_x>, <I_y> and <I_z>.";

  return TrussMe_FEM:-GetSpringStiffness([E*A/L, 0, 0], [G*(I_x+I_y)/L, 0, 0]);

end proc: # GetRodStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetBeamStiffness::static := proc(
  A::algebraic,
  E::algebraic,
  G::algebraic,
  L::algebraic,
  I_x::algebraic,
  I_y::algebraic,
  I_z::algebraic,
  $)::STIFFNESS;

  description "Get the stiffness matrix of a rod (only axial-stiffness) "
    "given the cross-section area <A>, the Young's modulus <E>, the shear "
    "modulus <G>, the length <L>, and the inertia of the cross-section "
    "<I_x>, <I_y> and <I_z>.";

  local EA_L, GJ_L, EIy_L, EIy_L2, EIy_L3, EIz_L, EIz_L2, EIz_L3, vec, mat;

  EA_L  := E*A/L;
  GJ_L  := G*(I_x+I_y)/L;
  EIy_L := E*I_y/L; EIy_L2 := E*I_y/L^2; EIy_L3 := E*I_y/L^3;
  EIz_L := E*I_z/L; EIz_L2 := E*I_z/L^2; EIz_L3 := E*I_z/L^3;

  return Matrix(
    <<EA_L, 0, 0, 0, 0, 0, -EA_L, 0, 0, 0, 0, 0>|
     <0, 12*EIz_L3, 0, 0, 0, 6*EIz_L2, 0, -12*EIz_L3, 0, 0, 0, 6*EIz_L2>|
     <0, 0, 12*EIy_L3, 0, -6*EIy_L2, 0, 0, 0, -12*EIy_L3, 0, -6*EIy_L2, 0>|
     <0, 0, 0, GJ_L, 0, 0, 0, 0, 0, -GJ_L, 0, 0>|
     <0, 0, -6*EIy_L2, 0, 4*EIy_L, 0, 0, 0, 6*EIy_L2, 0, 2*EIy_L, 0>|
     <0, 6*EIz_L2, 0, 0, 0, 4*EIz_L, 0, -6*EIz_L2, 0, 0, 0, 2*EIz_L>|
     <-EA_L, 0, 0, 0, 0, 0, EA_L, 0, 0, 0, 0, 0>|
     <0, -12*EIz_L3, 0, 0, 0, -6*EIz_L2, 0, 12*EIz_L3, 0, 0, 0, -6*EIz_L2>|
     <0, 0, -12*EIy_L3, 0, 6*EIy_L2, 0, 0, 0, 12*EIy_L3, 0, 6*EIy_L2, 0>|
     <0, 0, 0, -GJ_L, 0, 0, 0, 0, 0, GJ_L, 0, 0>|
     <0, 0, -6*EIy_L2, 0, 2*EIy_L, 0, 0, 0, 6*EIy_L2, 0, 4*EIy_L, 0>|
     <0, 6*EIz_L2, 0, 0, 0, 2*EIz_L, 0, -6*EIz_L2, 0, 0, 0, 4*EIz_L>>,
    storage = sparse);
end proc: # GetBeamStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsELEMENT::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of ELEMENT type.";

  return type(var, table) and evalb(var["type"] = ELEMENT);
end proc: # IsELEMENT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeElement::static := proc(
  name::string,
  RF::FRAME,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  K::STIFFNESS,{
  distance::algebraic := -1
  }, $)::ELEMENT;

  description "Make an element with name <name> on reference frame <RF>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with stiffness <K>. "
    "Optional nodes distance <distance> can be specified.";

  local node_1, node_2, dofs_1, dofs_2, L, Q;

  if evalb(nops(N1) = 1) and type(N1, NODE) then
    node_1 := N1;
    dofs_1 := [1, 1, 1, 1, 1, 1];
  elif evalb(nops(N1) = 2) and type(N1[1], NODE) and type(N1[2], DOFS) then
    node_1 := N1[1];
    dofs_1 := N1[2];
  else
    error("invalid node 1 detected.");
  end if;

  if evalb(nops(N2) = 1) and type(N2, NODE) then
    node_2 := N2;
    dofs_2 := [1, 1, 1, 1, 1, 1];
  elif evalb(nops(N2) = 2) and type(N2[1], NODE) and type(N2[2], DOFS) then
    node_2 := N2[1];
    dofs_2 := N2[2];
  else
    error("invalid node 2 detected.");
  end if;

  if evalb(distance <> -1) then
    L := distance;
  else
    L := TrussMe_FEM:-Norm2(node_1["coordinates"] - node_2["coordinates"]);
  end if;

  Q := Matrix(12, storage = sparse);
  R := TrussMe_FEM:-Rotation(RF);
  Q[1..3,   1..3  ] := R;
  Q[4..6,   4..6  ] := R;
  Q[7..9,   7..9  ] := R;
  Q[10..12, 10..12] := R;

  return table([
    "type"      = ELEMENT,
    "name"      = name,
    "id"        = TrussMe_FEM:-GenerateId(),
    "frame"     = RF,
    "node_1"    = node_1["id"],
    "dofs_1"    = dofs_1,
    "node_2"    = node_2["id"],
    "dofs_2"    = dofs_2,
    "length"    = L,
    "stiffness" = K,
    "rotation"  = Q
  ]);
end proc: # MakeElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeSpring::static := proc(
  name::string,
  RF::FRAME,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  {
  K::{algebraic, list(algebraic)} := 0,
  T::{algebraic, list(algebraic)} := 0,
  distance::algebraic             := -1
  }, $)::ELEMENT;

  description "Make a spring element with name <name> on reference frame <RF>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with translational "
    "stiffnesses <K> or <K_x, K_y, K_z>, and torsional stiffnesses <T> or "
    "<T_x, T_y, T_z>. Optional nodes distance <distance> can be specified.";

  return TrussMe_FEM:-MakeElement(name, RF, N1, N2,
    TrussMe_FEM:-GetSpringStiffness(K, T), parse("distance") = distance
  );
end proc: # MakeSpring

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeRod::static := proc(
  name::string,
  RF::FRAME,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  {
  material::MATERIAL       := TrussMe_FEM:-MakeCarbonSteel(),
  area::algebraic          := 0,
  inertia::list(algebraic) := [0, 0, 0],
  distance::algebraic      := -1
  }, $)::ELEMENT;

  description "Make a rod element with name <name> on reference frame <RF>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with material <material> and "
    "cross-section area <area> and inertia <inertia>. Optional nodes distance "
    "<distance> can be specified.";

  local L;

  if evalb(distance <> -1) then
    L := distance;
  else
    if type(N1, NODE) and type(N2, NODE) then
      L := TrussMe_FEM:-Norm2(N1["coordinates"] - N2["coordinates"]);
    elif evalb(nops(N1) = 2) and type(N1[1], NODE) and type(N1[2], DOFS) and
        evalb(nops(N2) = 2) and type(N2[1], NODE) and type(N2[2], DOFS) then
      L := TrussMe_FEM:-Norm2(N1[1]["coordinates"] - N2[2]["coordinates"]);
    else
      error("invalid nodes detected.");
    end if;
  end if;

  return TrussMe_FEM:-MakeElement(name, RF, N1, N2, TrussMe_FEM:-GetRodStiffness(
    area, material["elastic_modulus"], material["shear_modulus"], L, op(inertia)
  ), parse("distance") = distance);
end proc: # MakeRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBeam::static := proc(
  name::string,
  RF::FRAME,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  {
  material::MATERIAL       := TrussMe_FEM:-MakeCarbonSteel(),
  area::algebraic          := 0,
  inertia::list(algebraic) := [0, 0, 0],
  distance::algebraic      := -1
  }, $)::ELEMENT;

  description "Make a beam element with name <name> on reference frame <RF>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with material <material> and "
    "cross-section area <area> and inertia <inertia>. Optional nodes distance "
    "<distance> can be specified.";

  local L;

  if evalb(distance <> -1) then
    L := distance;
  else
    if type(N1, NODE) and type(N2, NODE) then
      L := TrussMe_FEM:-Norm2(N1["coordinates"] - N2["coordinates"]);
    elif evalb(nops(N1) = 2) and type(N1[1], NODE) and type(N1[2], DOFS) and
        evalb(nops(N2) = 2) and type(N2[1], NODE) and type(N2[2], DOFS) then
      L := TrussMe_FEM:-Norm2(N1[1]["coordinates"] - N2[2]["coordinates"]);
    else
      error("invalid nodes detected.");
    end if;
  end if;

  return TrussMe_FEM:-MakeElement(name, RF, N1, N2, TrussMe_FEM:-GetBeamStiffness(
    area, material["elastic_modulus"], material["shear_modulus"], L, op(inertia)
  ), parse("distance") = distance);
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSTRUCTURE::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of STRUCTURE type.";

  return type(var[1], list(NODE)) and type(var[2], list(ELEMENT));
end proc: # IsSTRUCTURE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNodesDofs::static := proc(
  node::list(NODE),
  $)::list;

  description "Get the list of dofs of the structure <s>.";

  local node;

  dofs := [];
  for node in s do
    if type(node, NODE) then
      dofs := [op(dofs), op(node["dofs"])];
    end if;
  end do;
  return dofs;
end proc: # GetDofs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export StiffnessTransformation::static := proc(
  nodes::list(NODE),
  $)::Matrix;

  description "Compute the transformation matrix of the nodes <nodes>.";

  local T, i;

  T := Matrix(12 * nops(nodes), shape = identity, storage = sparse);
  for i from 1 to nops(nodes) do
    if type(i, SUPPORT) then
      T[12*i-11..12*i-6, 12*i-11..12*i-6] := i["frame"];
      T[12*i-5..12*i,    12*i-5..12*i   ] := i["frame"];
    end if;
  end do;
  local T;
end proc: # Transformation

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GlobalStiffness::static := proc(
  nodes::list(NODE),
  elements::list(ELEMENT),
  $)::Matrix;

  description "Compute the global stiffness matrix of the nodes <nodes> and "
    "elements <elements>.";

  local K, i;

  K := Matrix(12 * nops(nodes), storage = sparse);
  for i from 1 to nops(elements) do
    j := TrussMe_FEM:-GetObjById(nodes, i["node_1"], parse("position") = true);
    k := TrussMe_FEM:-GetObjById(nodes, i["node_2"], parse("position") = true);
    dofs_1 := <op(i["dofs_1"]), op(i["dofs_1"])>;
    dofs_2 := <op(i["dofs_2"]), op(i["dofs_2"])>;
    K[12*j-11..12*j, 12*j-11..12*j] := i["stiffness"].dofs_1.i["rotation"];
    K[12*k-11..12*k, 12*k-11..12*k] := i["stiffness"].dofs_2.i["rotation"];
  end do;
  local K;
end proc: # GlobalStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GlobalStiffnessPrime::static := proc(
  nodes::list(NODE),
  elements::list(ELEMENT),
  $)::Matrix;

  description "Compute the global stiffness matrix of the nodes <nodes> and "
    "elements <elements>.";

  local T, K;

  T := TrussMe_FEM:-StiffnessTransformation(nodes);
  K := TrussMe_FEM:-GlobalStiffness(nodes, elements);
  return T.K.LinearAlgebra:-Transpose(T);
end proc: # GlobalStiffnessPrime

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

