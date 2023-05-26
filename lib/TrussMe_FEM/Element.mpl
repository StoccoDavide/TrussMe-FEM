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
  if type(K, list(algebraic)) and evalb(nops(K) = 3) then
    K_x := K[1]; K_y := K[2]; K_z := K[3];
  elif type(K, algebraic) then
    K_x := K; K_y := K; K_z := K;
  else
    error "<K> must be a list of 3 elements or a single element.";
  end if;

  # Retrieve the torsional stiffnesses
  if type(T, list(algebraic)) and evalb(nops(T) = 3) then
    T_x := T[1]; T_y := T[2]; T_z := T[3];
  elif type(T, algebraic) then
    T_x := T; T_y := T; T_z := T;
  else
    error "<T> must be a list of 3 elements or a single element.";
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
  I::list(algebraic),
  }, $)::STIFFNESS;

  description "Get the stiffness matrix of a rod (only axial-stiffness) "
    "given the cross-section area <A>, the Young's modulus <E>, the shear "
    "modulus <G>, the length <L>, and the inertia of the cross-section "
    "<I_x, I_y, I_z>.";

  if evalb(nops(I) <> 3) then
    error "<I> must be a list of 3 elements.";
  end if;

  return TrussMe_FEM:-GetSpringStiffness(
    parse("K_x") = E*A/L,           parse("K_y") = 0, parse("K_z") = 0,
    parse("T_x") = G*(I[1]+I[2])/L, parse("T_y") = 0, parse("T_z") = 0
  );

end proc: # GetRodStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetBeamStiffness::static := proc({
  A::algebraic,
  E::algebraic,
  G::algebraic,
  L::algebraic,
  I::list(algebraic),
  }, $)::STIFFNESS;

  description "Get the stiffness matrix of a rod (only axial-stiffness) "
    "given the cross-section area <A>, the Young's modulus <E>, the shear "
    "modulus <G>, the length <L>, and the inertia of the cross-section "
    "<I_x, I_y, I_z>.";

  local EA_L, GJ_L, EIy_L3, EIy_L2, EIy_L1, EIz_L3, EIz_L2, EIz_L1;

  if evalb(nops(I) <> 3) then
    error "<I> must be a list of 3 elements.";
  end if;

  EA_L := E*A/L; GJ_L := G*(I[1]+I[2])/L;
  EIy_L3 := E*I[2]/L^3; EIy_L2 := E*I[2]/L^2; EIy_L1 := E*I[2]/L;
  EIz_L3 := E*I[2]/L^3; EIz_L2 := E*I[2]/L^2; EIz_L1 := E*I[2]/L;
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

export MakeElement::static := proc(
  name::string,
  L::list({NODE, DOFS}),
  R::list({NODE, DOFS}),
  stiffness::STIFFNESS
  $)::ELEMENT;

  description "Make an element with name <name> and stiffness <stiffness>, "
    "connecting the dofs <L_dofs> on left node <L, L_dofs>, and connecting "
    "the dofs <R_dofs> on right node <R, R_dofs>.";

  if evalb(nops(L) <> 2) or not type(L[1], NODE) or not type(L[2], DOFS) then
    error "the first element of <L> must be a NODE and the second a DOFS";
  end if;
  if evalb(nops(R) <> 2) or not type(R[1], NODE) or not type(R[2], DOFS) then
    error "the first element of <R> must be a NODE and the second a DOFS";
  end if;

  # FIXME: lef node L may collide with length L change to something else

  return table([
    "type" = ELEMENT,
    "name" = name,
    "id"   = TrussMe_FEM:-GenerateId(),
    "L"    = [L[1]["id"], L[2]],
    "R"    = [R[1]["id"], R[2]],
    "K"    = stiffness
  ]);
end proc: # MakeElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeSpring::static := proc(
  name::string,
  L::list({NODE, DOFS}),
  R::list({NODE, DOFS}),
  {
  K::{algebraic, list(algebraic)},
  T::{algebraic, list(algebraic)}
  }, $)::ELEMENT;

  description "Make a spring element.";

  if evalb(nops(L) <> 2) or not type(L[1], NODE) or not type(L[2], DOFS) then
    error "the first element of <L> must be a NODE and the second a DOFS";
  end if;
  if evalb(nops(R) <> 2) or not type(R[1], NODE) or not type(R[2], DOFS) then
    error "the first element of <R> must be a NODE and the second a DOFS";
  end if;

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

export MakeRod::static := proc(
  name::string,
  L::list({NODE, DOFS}),
  R::list({NODE, DOFS}),
  {
  M::MATERIAL        := NULL,
  I::list(algebraic) := [0, 0, 0],
  A::algebraic       := 0,
  }, $)::ELEMENT;

  if evalb(nops(L) <> 2) or not type(L[1], NODE) or not type(L[2], DOFS) then
    error "the first element of <L> must be a NODE and the second a DOFS";
  end if;
  if evalb(nops(R) <> 2) or not type(R[1], NODE) or not type(R[2], DOFS) then
    error "the first element of <R> must be a NODE and the second a DOFS";
  end if;

  return table([
    "type" = ELEMENT,
    "name" = name,
    "id"   = TrussMe_FEM:-GenerateId(),
    "L"    = [L[1]["id"], L[2]],
    "R"    = [R[1]["id"], R[2]],
    "K"    = TrussMe_FEM:-GetRodStiffness(
      parse("M") = M, parse("I") = I, parse("A") = A
    );
  ]);

end proc: # MakeRod

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

