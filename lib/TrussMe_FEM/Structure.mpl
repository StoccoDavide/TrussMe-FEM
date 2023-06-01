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

export IsNODES::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of NODES type.";

  return type(var, list(NODE));
end proc: # IsNODES

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSUPPORT::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of SUPPORT type.";

  return type(var, NODE) and has(map(has, var["dofs"], 0), true);
end proc: # IsSUPPORT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSUPPORTS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of SUPPORTS type.";

  return type(var, list(SUPPORT));
end proc: # IsSUPPORTS

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
  coordinates::{list, POINT, VECTOR},
  {
  frame::FRAME                   := Matrix(4, shape = identity),
  dofs::DOFS                     := [1, 1, 1, 1, 1, 1],
  displacements::list(algebraic) := [0, 0, 0, 0, 0, 0]
  }, $)::NODE;

  description "Create a node with name <name> at coordinates <coordinates> in the "
    "reference frame <frame>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained to the ground "
    "in the direction given from <frame>.";

  local coordinates_tmp;

  if type(coordinates, list) and evalb(nops(coordinates) = 3) then
    coordinates_tmp := coordinates;
  elif type(coordinates, POINT) or type(coordinates, VECTOR) then
    coordinates_tmp := [TrussMe_FEM:-CompXYZ(coordinates)];
  else
    error("<coordinates> must be a list of 3 elements, a POINT or a VECTOR.");
  end if;

  if evalb(nops(displacements) <> 6) then
    error("<displacements> must be a list of 6 elements.");
  elif evalb(add(dofs *~ displacements) <> 0) then
    error("<displacements> must be defined only for constrained dofs.");
  end if;

  return table([
    "type"                 = NODE,
    "name"                 = name,
    "id"                   = TrussMe_FEM:-GenerateId(),
    "frame"                = frame,
    "coordinates"          = coordinates_tmp,
    "dofs"                 = dofs,
    "displacements"        = displacements,
    "output_reactions"     = [], # Output reactions in the node frame
    "output_displacements" = []  # Output displacements in the node frame
  ]);
end proc: # MakeNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCompliantNode::static := proc(
  name::string,
  coordinates::{list, POINT, VECTOR},
  {
  frame::FRAME                    := Matrix(4, shape = identity),
  dofs::DOFS                      := [1, 1, 1, 1, 1, 1],
  K::{algebraic, list(algebraic)} := 0,
  T::{algebraic, list(algebraic)} := 0
  }, $)::NODE, ELEMENT, NODE;

  description "Create a node with name <name> at coordinates <coordinates> in the "
    "reference frame <frame>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained to the ground "
    "in the direction given from <frame>. The node is also connected to a "
    "compliant spring element with traslational stiffness <K> and torsional "
    "stiffness <T>.";

  local fixed, spring, compliant;

  fixed := TrussMe_FEM:-MakeNode(
    name, coordinates, parse("frame") = frame, parse("dofs") = dofs
  );
  compliant := TrussMe_FEM:-MakeNode(
    cat(name, "_compliant"), coordinates, parse("frame") = frame
  );
  spring    := TrussMe_FEM:-MakeSpring(
    cat(name, "_spring"), fixed, compliant,
    parse("frame") = frame, parse("K") = K, parse("T") = T
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
    error("<K> must be a list of 3 elements or a single element.");
  end if;

  if type(T, list(algebraic)) and evalb(nops(T) = 3) then
    T_x := T[1]; T_y := T[2]; T_z := T[3];
  elif type(T, algebraic) then
    T_x := T; T_y := T; T_z := T;
  else
    error("<T> must be a list of 3 elements or a single element.");
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

export IsELEMENTS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of ELEMENTS type.";

  return type(var, list(ELEMENT));
end proc: # IsELEMENTS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeElement::static := proc(
  name::string,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  K::STIFFNESS,
  {
  frame::FRAME        := TrussMe_FEM:-GenerateGenericFrame(name),
  distance::algebraic := -1
  }, $)::ELEMENT;

  description "Make an element with name <name> on reference frame <frame>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with stiffness <K>. "
    "Optional nodes distance <distance> can be specified.";

  local node_1, node_2, dofs_1, dofs_2, L, R, Q;

  if evalb(nops(N1) = 1) and type(N1, NODE) then
    node_1 := N1;
    dofs_1 := [0, 0, 0, 0, 0, 0];
  elif evalb(nops(N1) = 2) and type(N1[1], NODE) and type(N1[2], DOFS) then
    node_1 := N1[1];
    dofs_1 := N1[2];
  else
    error("invalid node 1 detected.");
  end if;

  if evalb(nops(N2) = 1) and type(N2, NODE) then
    node_2 := N2;
    dofs_2 := [0, 0, 0, 0, 0, 0];
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
  R := TrussMe_FEM:-Rotation(frame);
  Q[1..3,   1..3  ] := R;
  Q[4..6,   4..6  ] := R;
  Q[7..9,   7..9  ] := R;
  Q[10..12, 10..12] := R;

  return table([
    "type"      = ELEMENT,
    "name"      = name,
    "id"        = TrussMe_FEM:-GenerateId(),
    "frame"     = frame,
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
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  {
  K::{algebraic, list(algebraic)} := 0,
  T::{algebraic, list(algebraic)} := 0,
  frame::FRAME                    := TrussMe_FEM:-GenerateGenericFrame(name),
  distance::algebraic             := -1
  }, $)::ELEMENT;

  description "Make a spring element with name <name> on reference frame <frame>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with translational "
    "stiffnesses <K> or <K_x, K_y, K_z>, and torsional stiffnesses <T> or "
    "<T_x, T_y, T_z>. Optional nodes distance <distance> can be specified.";

  return TrussMe_FEM:-MakeElement(
    name, N1, N2, TrussMe_FEM:-GetSpringStiffness(K, T),
    parse("frame") = frame, parse("distance") = distance
  );
end proc: # MakeSpring

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeRod::static := proc(
  name::string,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  {
  material::MATERIAL       := TrussMe_FEM:-MakeCarbonSteel(),
  area::algebraic          := 0,
  inertia::list(algebraic) := [0, 0, 0],
  frame::FRAME             := TrussMe_FEM:-GenerateGenericFrame(name),
  distance::algebraic      := -1
  }, $)::ELEMENT;

  description "Make a rod element with name <name> on reference frame <frame>, "
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

  return TrussMe_FEM:-MakeElement(name, N1, N2, TrussMe_FEM:-GetRodStiffness(
    area, material["elastic_modulus"], material["shear_modulus"], L, op(inertia)
  ), parse("frame") = frame, parse("distance") = distance);
end proc: # MakeRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBeam::static := proc(
  name::string,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  {
  material::MATERIAL       := TrussMe_FEM:-MakeCarbonSteel(),
  area::algebraic          := 0,
  inertia::list(algebraic) := [0, 0, 0],
  frame::FRAME             := TrussMe_FEM:-GenerateGenericFrame(name),
  distance::algebraic      := -1
  }, $)::ELEMENT;

  description "Make a beam element with name <name> on reference frame <frame>, "
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

  return TrussMe_FEM:-MakeElement(name, N1, N2, TrussMe_FEM:-GetBeamStiffness(
    area, material["elastic_modulus"], material["shear_modulus"], L, op(inertia)
  ), parse("frame") = frame, parse("distance") = distance);
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSTRUCTURE::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of STRUCTURE type.";

  return type(var[1], NODES) and type(var[2], ELEMENTS);
end proc: # IsSTRUCTURE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNodalDofs::static := proc(
  nodes::NODES,
  $)::Vector;

  description "Get the vector of nodal dofs of nodes <nodes>.";

  local dofs, i;

  printf("TrussMe_FEM:-GetNodalDofs() ... ");
  dofs := [];
  for i in nodes do
    if type(i, NODE) then
      dofs := [op(dofs), op(i["dofs"])];
    end if;
  end do;
  printf("DONE\n");
  return convert(dofs, Vector);
end proc: # GetNodalDofs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNodalDisplacements::static := proc(
  nodes::NODES,
  $)::Vector;

  description "Get the vector of nodal displacements of nodes <nodes>.";

  local displacements, i;

  displacements := [];
  for i in nodes do
    if type(i, NODE) then
      displacements := [op(displacements), op(i["displacements"])];
    end if;
  end do;
  return convert(displacements, Vector);
end proc: # GetNodalDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export StiffnessTransformation::static := proc(
  nodes::NODES,
  $)::Matrix;

  description "Compute the transformation matrix of the nodes <nodes>.";

  local T, i;

  T := Matrix(6*nops(nodes), (i,j) -> `if`(evalb(i=j), 1, 0), storage = sparse);
  for i from 1 to nops(nodes) do
    T[6*i-5..6*i-3, 6*i-5..6*i-3] := nodes[i]["frame"][1..3, 1..3];
    T[6*i-2..6*i,   6*i-2..6*i]   := nodes[i]["frame"][1..3, 1..3];
  end do;
  return T;
end proc: # StiffnessTransformation

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GlobalStiffness::static := proc(
  nodes::NODES,
  elements::ELEMENTS,
  $)::Matrix;

  description "Compute the global stiffness matrix of the nodes <nodes> and "
    "elements <elements>.";

  local K, element_k, i, j, k;

  printf("TrussMe_FEM:-GlobalStiffness() ... ");
  K := Matrix(6*nops(nodes), storage = sparse);
  for i from 1 to nops(elements) do
    # Nodes positions
    j := TrussMe_FEM:-GetObjById(nodes, elements[i]["node_1"], parse("position") = true);
    k := TrussMe_FEM:-GetObjById(nodes, elements[i]["node_2"], parse("position") = true);
    # Element stiffness contribution selecting only constrained dofs (= 0)
    element_k := elements[i]["stiffness"].LinearAlgebra:-DiagonalMatrix(
      <op(-(elements[i]["dofs_1"] -~ 1)), op(-(elements[i]["dofs_2"] -~ 1))>
    );
    element_k := elements[i]["rotation"].element_k;
    K[6*j-5..6*j, 6*j-5..6*j] := K[6*j-5..6*j, 6*j-5..6*j] + element_k[1..6,  1..6 ];
    K[6*j-5..6*j, 6*k-5..6*k] := K[6*j-5..6*j, 6*k-5..6*k] + element_k[1..6,  7..12];
    K[6*k-5..6*k, 6*j-5..6*j] := K[6*k-5..6*k, 6*j-5..6*j] + element_k[7..12, 1..6 ];
    K[6*k-5..6*k, 6*k-5..6*k] := K[6*k-5..6*k, 6*k-5..6*k] + element_k[7..12, 7..12];
  end do;
  printf("DONE\n");
  return K;
end proc: # GlobalStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GlobalStiffnessPrime::static := proc(
  nodes::NODES,
  elements::ELEMENTS,
  $)::Matrix;

  description "Compute the global stiffness matrix of the nodes <nodes> and "
    "elements <elements>.";

  local T, K;

  printf("TrussMe_FEM:-GlobalStiffnessPrime() ... ");
  T := TrussMe_FEM:-StiffnessTransformation(nodes);
  K := TrussMe_FEM:-GlobalStiffness(nodes, elements);
  printf("DONE\n");
  return T.K.LinearAlgebra:-Transpose(T);
end proc: # GlobalStiffnessPrime

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsFEM::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of FEM type.";

  return type(var, table) and evalb(var["type"] = FEM);
end proc: # IsFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateFEM::static := proc(
  nodes::NODES,
  elements::ELEMENTS,
  loads::LOADS,
  $)::FEM;

  description "Generate a FEM structure from the nodes <nodes>, elements "
    "<elements>, and loads <loads>.";

  return table([
    "type"                 = FEM,
    "nodes"                = map(x -> x["id"], nodes),
    "elements"             = map(x -> x["id"], elements),
    "dofs"                 = TrussMe_FEM:-GetNodalDofs(nodes),
    "displacements"        = TrussMe_FEM:-GetNodalDisplacements(nodes),
    "stiffness"            = TrussMe_FEM:-GlobalStiffnessPrime(nodes, elements),
    "loads"                = TrussMe_FEM:-GetNodalLoads(nodes, loads),
    "output_reactions"     = [], # Output reactions in the ground frame
    "output_displacements" = []  # Output displacements in the ground frame
  ]);
end proc: # GenerateFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CheckNodes::static := proc(
  fem::FEM,
  nodes::NODES,
  $)::boolean;

  description "Check if the nodes <nodes> associated to the FEM structure "
    "<fem> are valid (subset of fem['nodes']).";

  local fem_nodes, set_nodes;

  fem_nodes := convert(fem["nodes"], set);
  set_nodes := convert(map(x -> x["id"], nodes), set);

  return evalb(nops(set_nodes union fem_nodes) = nops(fem_nodes));
end proc: # CheckNodes

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CheckElements::static := proc(
  fem::FEM,
  elements::ELEMENTS,
  $)::boolean;

  description "Check if the elements <elements> associated to the FEM structure "
    "<fem> are valid (subset of fem['elements']).";

  local fem_elements, set_elements;

  fem_elements := convert(fem["elements"], set);
  set_elements := convert(map(x -> x["id"], elements), set);

  return evalb(nops(set_elements union fem_elements) = nops(fem_elements));
end proc: # CheckElements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SplitFEM::static := proc(
  fem::FEM,
  $)

  local d, K, f, perm, n, K_ff, K_fs, K_sf, K_ss, d_s, f_f;

  # Apply permutation to stiffness and loads to split free and constrained dofs
  perm := sort(fem["dofs"], `>`, output = 'permutation');
  K := fem["stiffness"][perm, perm];
  d := fem["displacements"][perm];
  f := fem["loads"][perm];

  # Compute output displacements and reactions
  n := add(fem["dofs"]);
  return table([
    "perm" = perm,
    "K_ff" = K[1..n, 1..n],
    "K_fs" = K[1..n, n+1..-1],
    "K_sf" = K[n+1..-1, 1..n],
    "K_ss" = K[n+1..-1, n+1..-1],
    "d_s"  = d[n+1..-1],
    "f_f"  = f[1..n]
  ]);
end proc: # SplitFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SolveFEM::static := proc(
  fem::FEM,
  $)

  description "Solve the FEM structure <fem>.";

  local split_FEM, i, mrep;

  # Split free and constrained dofs
  split_FEM := TrussMe_FEM:-SplitFEM(fem);

  # Solve free dofs
  split_FEM["d_f"] := LinearAlgebra:-LinearSolve(
    split_FEM["K_ff"], split_FEM["f_f"] - split_FEM["K_fs"].split_FEM["d_s"]
  );

  # Solve reactions
  split_FEM["f_s"] :=
    split_FEM["K_sf"].split_FEM["d_f"] - split_FEM["K_ss"].split_FEM["d_s"];

  # Store output displacements and reactions and restore initial permutation
  mrep := [seq(i, i = 1..nops(split_FEM["perm"]))];
  for i from 1 to nops(mrep) do
    mrep[split_FEM["perm"][i]] := i;
  end do;
  fem["output_displacements"] := <split_FEM["d_f"], split_FEM["d_s"]>[mrep];
  fem["output_reactions"]     := <split_FEM["f_f"], split_FEM["f_s"]>[mrep];

  if evalb(_nresults = 1) then
    return split_FEM;
  else
    return NULL;
  end if;
end proc: # SolveFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export StoreFEM::static := proc(
  fem::FEM,
  nodes::NODES,
  $)

  description "Store the FEM structure <fem> in the nodes <nodes>.";

  local i, tmp;

  if not TrussMe_FEM:-CheckNodes(fem, nodes) then
    error("invalid nodes detected.");
  end if;

  for i from 1 to nops(nodes) do
    tmp := nodes[i];
    tmp["output_displacements"] := fem["output_displacements"][6*i-5..6*i];
    tmp["output_reactions"]     := fem["output_reactions"][6*i-5..6*i];
  end do;
  return NULL;
end proc: # StoreFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
