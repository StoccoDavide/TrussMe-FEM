# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                    _____ _                           _                      #
#                   | ____| | ___ _ __ ___   ___ _ __ | |_                    #
#                   |  _| | |/ _ \ '_ ` _ \ / _ \ '_ \| __|                   #
#                   | |___| |  __/ | | | | |  __/ | | | |_                    #
#                   |_____|_|\___|_| |_| |_|\___|_| |_|\__|                   #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSTIFFNESS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of STIFFNESS type.";

  return type(var, Matrix) and evalb(var["type"] = STIFFNESS) and
    evalb(LinearAlgebra:-RowDimension(var) = 12) and
    evalb(LinearAlgebra:-ColumnDimension(var) = 12);
end proc: # IsSTIFFNESS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetSpringStiffness::static := proc({
  K::{algebraic, list(algebraic)},
  T::{algebraic, list(algebraic)},
  }, $)::STIFFNESS;

  description "Get the stiffness matrix of a spring given the translational "
    "stiffnesses <K> or <K_x, K_y, K_z>, and the torsional stiffnesses <T> or "
    "<T_x, T_y, T_z>.";

  local K_x, K_y, K_z, T_x, T_y, T_z;

  # Retrieve the translational stiffnesses
  if type(K, list(algebraic)) then
    K_x := K[1]; K_y := K[2]; K_z := K[3];
  else
    K_x := K; K_y := K; K_z := K;
  end if;

  # Retrieve the torsional stiffnesses
  if type(T, list(algebraic)) then
    T_x := T[1]; T_y := T[2]; T_z := T[3];
  else
    T_x := T; T_y := T; T_z := T;
  end if;

  # Return the stiffness matrix
  return <
    <K_x, 0, 0, 0, 0, 0, -K_x, 0, 0, 0, 0, 0>,
    <0, K_y, 0, 0, 0, 0, 0, -K_y, 0, 0, 0, 0>,
    <0, 0, K_z, 0, 0, 0, 0, 0, -K_z, 0, 0, 0>,
    <0, 0, 0, T_x, 0, 0, 0, 0, 0, -T_x, 0, 0>,
    <0, 0, 0, 0, T_y, 0, 0, 0, 0, 0, -T_y, 0>,
    <0, 0, 0, 0, 0, T_z, 0, 0, 0, 0, 0, -T_z>,
    <-K_x, 0, 0, 0, 0, 0, K_x, 0, 0, 0, 0, 0>,
    <0, -K_y, 0, 0, 0, 0, 0, K_y, 0, 0, 0, 0>,
    <0, 0, -K_z, 0, 0, 0, 0, 0, K_z, 0, 0, 0>,
    <0, 0, 0, -T_x, 0, 0, 0, 0, 0, T_x, 0, 0>,
    <0, 0, 0, 0, -T_y, 0, 0, 0, 0, 0, T_y, 0>,
    <0, 0, 0, 0, 0, -T_z, 0, 0, 0, 0, 0, T_z>
  >;
end proc: # GetSpringStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetRodStiffness::static := proc({
  A::algebraic,
  E::algebraic,
  G::algebraic,
  L::algebraic,
  J::algebraic,
  }, $)::STIFFNESS;

  description "Get the 12x12 stiffness matrix of a rod (only axial-stiffness) "
    "given the cross-section area <A>, the Young's modulus <E>, the shear "
    "modulus <G>, the length <L>, and the torsional constant <J>.";

  local EA_L, GJ_L, EIy_L3, EIy_L2, EIy_L1, EIz_L3, EIz_L2, EIz_L1;

  EA_L = E*A/L; GJ_L = G*J/L;
  return TrussMe_FEM:-GetSpringStiffness(
    parse("K_x") = EA_L, parse("K_y") = 0, parse("K_z") = 0,
    parse("T_x") = GJ_L, parse("T_y") = 0, parse("T_z") = 0
  );

end proc: # GetRodStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetBeamStiffness::static := proc({
  A::algebraic,
  E::algebraic,
  G::algebraic,
  L::algebraic,
  I_y::list(algebraic),
  I_z::list(algebraic),
  J::algebraic,
  }, $)::STIFFNESS;

  description "Get the stiffness matrix o";

  local EA_L, GJ_L, EIy_L3, EIy_L2, EIy_L1, EIz_L3, EIz_L2, EIz_L1;

  EA_L = E*A/L; GJ_L = G*J/L;
  EIy_L3 := E*I_y/L^3; EIy_L2 := E*I_y/L^2; EIy_L1 := E*I_y/L;
  EIz_L3 := E*I_z/L^3; EIz_L2 := E*I_z/L^2; EIz_L1 := E*I_z/L;
  return <
    <EA_L, 0, 0, 0, 0, 0, -EA_L, 0, 0, 0, 0, 0>,
    <0, 12*EIz_L3, 0, 0, 0, 6*EIz_L2, 0, -12*EIz_L3, 0, 0, 0, 6*EIz_L2>,
    <0, 0, 12*EIy_L3, 0, -6*EIy_L2, 0, 0, 0, -12*EIy_L3, 0, -6*EIy_L2, 0>,
    <0, 0, 0, GJ_L, 0, 0, 0, 0, 0, -GJ_L, 0, 0>,
    <0, 0, -6*EIy_L2, 0, 4*EIy_L1, 0, 0, 0, 6*EIy_L2, 0, 2*EIy_L1, 0>,
    <0, 6*EIz_L2, 0, 0, 0, 4*EIz_L1, 0, -6*EIz_L2, 0, 0, 0, 2*EIz_L1>,
    <-EA_L, 0, 0, 0, 0, 0, EA_L, 0, 0, 0, 0, 0>,
    <0, -12*EIz_L3, 0, 0, 0, -6*EIz_L2, 0, 12*EIz_L3, 0, 0, 0, -6*EIz_L2>,
    <0, 0, -12*EIy_L3, 0, 6*EIy_L2, 0, 0, 0, 12*EIy_L3, 0, 6*EIy_L2, 0>,
    <0, 0, 0, -GJ_L, 0, 0, 0, 0, 0, GJ_L, 0, 0>,
    <0, 0, -6*EIy_L2, 0, 2*EIy_L1, 0, 0, 0, 6*EIy_L2, 0, 4*EIy_L1, 0>,
    <0, 6*EIz_L2, 0, 0, 0, 2*EIz_L1, 0, -6*EIz_L2, 0, 0, 0, 4*EIz_L1>
  >;
end proc: # GetBeamStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsELEMENT::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of ELEMENT type.";

  return type(var, table) and evalb(var["type"] = ELEMENT);
end proc: # IsELEMENT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeElement::static := proc({
  name::string,
  node_L::NODE,
  dofs_L::list(boolean),
  node_R::NODE,
  dofs_R::list(boolean),
  material::MATERIAL,
  stiffness::STIFFNESS
  $)::ELEMENT;

  return table([
    "type"      = ELEMENT,
    "name"      = name,
    "id"        = TrussMe_FEM:-GenerateId(),
    "node_L"    = [node_L["id"], dofs_L],
    "node_R"    = [node_R["id"], dofs_R],
    "stiffness" = stiffness
  ]);
end proc: # MakeElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeSpring::static := proc(
  name::string,
  L::{NODE, list(integer)},
  R::{NODE, list(integer)},
  {
  K::{algebraic, list(algebraic)},
  T::{algebraic, list(algebraic)}
  }, $)::ELEMENT;

  return table([
    "type" = ELEMENT,
    "name" = name,
    "id"   = TrussMe_FEM:-GenerateId(),
    "L"    = [L[1]["id"], L[2]],
    "R"    = [R[1]["id"], R[2]],
    "K"    = TrussMe_FEM:-GetSpringStiffness(parse("K") = K, parse("T") = T);
  ]);
end proc: # MakeSpring

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeRod::static := proc({
  name::string,
  node_L::NODE,
  dofs_L::list(boolean),
  node_R::NODE,
  dofs_R::list(boolean),
  material::MATERIAL := NULL,
  inertia::INERTIA   := [0, 0, 0],
  inertia::INERTIA   := [0, 0, 0],
  }, $)::ELEMENT;

end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBeam::static := proc({
  name::string,
  node_L::NODE,
  dofs_L::list(boolean),
  node_R::NODE,
  dofs_R::list(boolean),
  material::MATERIAL := NULL,
  inertia::INERTIA   := [0, 0, 0],
  inertia::INERTIA   := [0, 0, 0],
  }, $)::ELEMENT;


  local node_1_tmp, node_2_tmp;

  return table([
    "type"   = ELEMENT,
    "name"   = name,
    "id"     = TrussMe_FEM:-GenerateId(),
    "node_L" = [node_L["id"], dofs_L],
    "node_R" = [node_R["id"], dofs_R],
    "material" = material
    "stiffness" = NULL,
  ]);
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

