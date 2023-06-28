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

export IsNODE := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of NODE type.";

  return type(var, table) and evalb(var["type"] = NODE);
end proc: # IsNODE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsNODES := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of NODES type.";

  return type(var, list(NODE));
end proc: # IsNODES

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSUPPORT := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of SUPPORT type.";

  return type(var, NODE) and has(map(has, var["dofs"], 0), true);
end proc: # IsSUPPORT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSUPPORTS := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of SUPPORTS type.";

  return type(var, list(SUPPORT));
end proc: # IsSUPPORTS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsDOFS := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of DOFS type.";

  return type(var, Vector) and evalb(LinearAlgebra:-Dimension(var) = 6) and
    not has(map(has, var, {0, 1}), false);
end proc: # IsDOFS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeNode := proc(
  name::string,
  coordinates::{Vector, POINT, VECTOR},
  {
  frame::FRAME                     := Matrix(4, shape = identity),
  dofs:: DOFS                      := <1, 1, 1, 1, 1, 1>,
  displacements::Vector(algebraic) := <0, 0, 0, 0, 0, 0>
  }, $)::NODE;

  description "Create a node with name <name> at coordinates <coordinates> in the "
    "reference frame <frame>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained in the "
    "direction given from (global to local) <frame>.";

  local coordinates_tmp, displacements_tmp;

  if type(coordinates, Vector) and
    evalb(LinearAlgebra:-Dimension(coordinates) = 3) then
    coordinates_tmp := coordinates;
  elif type(coordinates, POINT) or type(coordinates, VECTOR) then
    coordinates_tmp := coordinates[1..3];
  else
    error("<coordinates> must be a vector of 3 or 4 elements.");
  end if;

  if type(displacements, Vector) and
    evalb(LinearAlgebra:-Dimension(displacements) = 6) then
    displacements_tmp := displacements;
  else
    error("<displacements> must be a vector of 6 or 4 elements.");
  end if;

  if evalb(add(dofs *~ displacements_tmp) <> 0) then
    error("<displacements> must be defined only for constrained dofs.");
  end if;

  return table([
    "type"                 = NODE,
    "name"                 = name,
    "id"                   = TrussMe:-FEM:-GenerateId(),
    "frame"                = frame, # Reference frame from global to local
    "coordinates"          = coordinates_tmp, # Coordinates in the global frame
    "dofs"                 = dofs, # Constrained dofs in the local frame
    "displacements"        = displacements_tmp, # Displacements in the local frame
    "output_reactions"     = [], # Output reactions in the local frame
    "output_displacements" = []  # Output displacements in the local frame
  ]);
end proc: # MakeNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCompliantNode := proc(
  name::string,
  coordinates::{Vector, POINT, VECTOR},
  {
  frame::FRAME                     := Matrix(4, shape = identity),
  dofs:: DOFS                      := <1, 1, 1, 1, 1, 1>,
  displacements::Vector(algebraic) := <0, 0, 0, 0, 0, 0>,
  K::{algebraic, list(algebraic)}  := 0,
  T::{algebraic, list(algebraic)}  := 0
  }, $)::NODE, ELEMENT, NODE;

  description "Create a node with name <name> at coordinates <coordinates> in the "
    "reference frame <frame>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained in the "
    "direction given from (global to local) <frame>. The node is also connected "
    "to a compliant spring element with traslational stiffness <K> and torsional "
    "stiffness <T>.";

  local fixed, spring, compliant;

  fixed := TrussMe:-FEM:-MakeNode(
    name, coordinates, parse("frame") = frame, parse("dofs") = dofs,
    parse("displacements") = displacements
  );
  compliant := TrussMe:-FEM:-MakeNode(
    cat(name, "_compliant"), coordinates, parse("frame") = frame
  );
  spring := TrussMe:-FEM:-MakeSpring(
    cat(name, "_spring"), fixed, compliant,
    parse("frame") = frame, parse("K") = K, parse("T") = T
  );
  return fixed, spring, compliant;
end proc: # MakeCompliantNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSTIFFNESS := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of STIFFNESS type.";

  return type(var, Matrix) and
    evalb(LinearAlgebra:-RowDimension(var) = LinearAlgebra:-ColumnDimension(var)) and
    evalb(irem(LinearAlgebra:-RowDimension(var), 6) = 0) and
    evalb(irem(LinearAlgebra:-ColumnDimension(var), 6) = 0);
end proc: # IsSTIFFNESS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetSpringStiffness := proc(
  K::{algebraic, list(algebraic)},
  T::{algebraic, list(algebraic)},
  $)::STIFFNESS;

  description "Get the stiffness matrix of a spring given the translational "
    "stiffnesses <K> or <K_x, K_y, K_z>, and the torsional stiffnesses <T> or "
    "<T_x, T_y, T_z>.";

  local out, K_x, K_y, K_z, T_x, T_y, T_z;

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

  if TrussMe:-FEM:-m_WarningMode then
    if type(K_x, numeric) and evalb(K_x < 0) then
      WARNING("negative x-axis translational stiffness detected.");
    end if;
    if type(K_y, numeric) and evalb(K_y < 0) then
      WARNING("negative y-axis translational stiffness detected.");
    end if;
    if type(K_z, numeric) and evalb(K_z < 0) then
      WARNING("negative z-axis translational stiffness detected.");
    end if;
    if type(T_x, numeric) and evalb(T_x < 0) then
      WARNING("negative x-axis torsional stiffness detected.");
    end if;
    if type(T_y, numeric) and evalb(T_y < 0) then
      WARNING("negative y-axis torsional stiffness detected.");
    end if;
    if type(T_z, numeric) and evalb(T_z < 0) then
      WARNING("negative z-axis torsional stiffness detected.");
    end if;
  end if;

  out := Matrix(12, storage = sparse);
  out[1, 1]   :=  K_x; out[2, 2]   :=  K_y; out[3, 3]   :=  K_z;
  out[4, 4]   :=  T_x; out[5, 5]   :=  T_y; out[6, 6]   :=  T_z;
  out[7, 7]   :=  K_x; out[8, 8]   :=  K_y; out[9, 9]   :=  K_z;
  out[10, 10] :=  T_x; out[11, 11] :=  T_y; out[12, 12] :=  T_z;
  out[1, 7]   := -K_x; out[2, 8]   := -K_y; out[3, 9]   := -K_z;
  out[7, 1]   := -K_x; out[8, 2]   := -K_y; out[9, 3]   := -K_z;
  out[4, 10]  := -T_x; out[5, 11]  := -T_y; out[6, 12]  := -T_z;
  out[10, 4]  := -T_x; out[11, 5]  := -T_y; out[12, 6]  := -T_z;
  return out;
end proc: # GetSpringStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetRodStiffness := proc(
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
    "modulus <G>, the length <L>, and the cross-section inertia <I_x>, <I_y> "
    "and <I_z>.";

  if TrussMe:-FEM:-m_WarningMode then
    if type(A, numeric) and evalb(A <= 0) then
      WARNING("null or negative cross-section area detected.");
    end if;
    if type(E, numeric) and evalb(E <= 0) then
      WARNING("null or negative elastic modulus detected.");
    end if;
    if type(G, numeric) and evalb(G <= 0) then
      WARNING("null or negative shear modulus detected.");
    end if;
    if type(L, numeric) and evalb(L <= 0) then
      WARNING("null or negative length detected.");
    end if;
    if evalb(I_x <> 0) then
      WARNING("non null x-axis inertia detected.");
    end if;
    if type(I_y, numeric) and evalb(I_y <= 0) then
      WARNING("null or negative y-axis inertia detected.");
    end if;
    if type(I_z, numeric) and evalb(I_z <= 0) then
      WARNING("null or negative z-axis inertia detected.");
    end if;
  end if;

  return TrussMe:-FEM:-GetSpringStiffness([E*A/L, 0, 0], [G*(I_y+I_z)/L, 0, 0]);

end proc: # GetRodStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetBeamStiffness := proc(
  A::algebraic,
  E::algebraic,
  G::algebraic,
  L::algebraic,
  I_x::algebraic,
  I_y::algebraic,
  I_z::algebraic,
  $)::STIFFNESS;

  description "Get the stiffness matrix of a lean beam given the cross-section "
    "area <A>, the Young's modulus <E>, the shear modulus <G>, the length <L>, "
    "and the inertia of the cross-section <I_x>, <I_y> and <I_z>.";

  local out, EA_L, GJ_L, EIy_L, EIy_L2, EIy_L3, EIz_L, EIz_L2, EIz_L3;

  if TrussMe:-FEM:-m_WarningMode then
    if type(A, numeric) and evalb(A <= 0) then
      WARNING("null or negative cross-section area detected.");
    end if;
    if type(E, numeric) and evalb(E <= 0) then
      WARNING("null or negative elastic modulus detected.");
    end if;
    if type(G, numeric) and evalb(G <= 0) then
      WARNING("null or negative shear modulus detected.");
    end if;
    if type(L, numeric) and evalb(L <= 0) then
      WARNING("null or negative length detected.");
    end if;
    if type(I_x, numeric) and evalb(I_x <= 0) then
      WARNING("null or negative x-axis inertia detected.");
    end if;
    if type(I_y, numeric) and evalb(I_y <= 0) then
      WARNING("null or negative y-axis inertia detected.");
    end if;
    if type(I_z, numeric) and evalb(I_z <= 0) then
      WARNING("null or negative z-axis inertia detected.");
    end if;
  end if;

  EA_L  := E*A/L;
  GJ_L  := G*(I_y+I_z)/L;
  EIy_L := E*I_y/L; EIy_L2 := EIy_L/L; EIy_L3 := EIy_L2/L;
  EIz_L := E*I_z/L; EIz_L2 := EIz_L/L; EIz_L3 := EIz_L2/L;

  out := Matrix(12, storage = sparse);
  out[1, 1]  := EA_L;       out[1, 7]   := -EA_L;
  out[2, 2]  := 12*EIz_L3;  out[2, 6]   := 6*EIz_L2;
  out[2, 8]  := -12*EIz_L3; out[2, 12]  := 6*EIz_L2;
  out[3, 3]  := 12*EIy_L3;  out[3, 5]   := -6*EIy_L2;
  out[3, 9]  := -12*EIy_L3; out[3, 11]  := -6*EIy_L2;
  out[4, 4]  := GJ_L;       out[4, 10]  := -GJ_L;
  out[5, 3]  := -6*EIy_L2;  out[5, 5]   := 4*EIy_L;
  out[5, 9]  := 6*EIy_L2;   out[5, 11]  := 2*EIy_L;
  out[6, 2]  := 6*EIz_L2;   out[6, 6]   := 4*EIz_L;
  out[6, 8]  := -6*EIz_L2;  out[6, 12]  := 2*EIz_L;
  out[7, 1]  := -EA_L;      out[7, 7]   := EA_L;
  out[8, 2]  := -12*EIz_L3; out[8, 6]   := -6*EIz_L2;
  out[8, 8]  := 12*EIz_L3;  out[8, 12]  := -6*EIz_L2;
  out[9, 3]  := -12*EIy_L3; out[9, 5]   := 6*EIy_L2;
  out[9, 9]  := 12*EIy_L3;  out[9, 11]  := 6*EIy_L2;
  out[10, 4] := -GJ_L;      out[10, 10] := GJ_L;
  out[11, 3] := -6*EIy_L2;  out[11, 5]  := 2*EIy_L;
  out[11, 9] := 6*EIy_L2;   out[11, 11] := 4*EIy_L;
  out[12, 2] := 6*EIz_L2;   out[12, 6]  := 2*EIz_L;
  out[12, 8] := -6*EIz_L2;  out[12, 12] := 4*EIz_L;
  return out;
end proc: # GetBeamStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetTimoshenkoBeamStiffness := proc(
  A::algebraic,
  E::algebraic,
  G::algebraic,
  L::algebraic,
  I_x::algebraic,
  I_y::algebraic,
  I_z::algebraic,
  $)::STIFFNESS;

  description "Get the stiffness matrix of a Timoshenko's (thick) beam "
    "given the cross-section area <A>, the Young's modulus <E>, the shear "
    "modulus <G>, the length <L>, and the inertia of the cross-section "
    "<I_x>, <I_y> and <I_z>.";

  local out, P_y, P_z, EA_L, GJ_L, EIy_L, EIy_L2, EIy_L3, EIz_L, EIz_L2, EIz_L3;

  if TrussMe:-FEM:-m_WarningMode then
    if type(A, numeric) and evalb(A <= 0) then
      WARNING("null or negative cross-section area detected.");
    end if;
    if type(E, numeric) and evalb(E <= 0) then
      WARNING("null or negative elastic modulus detected.");
    end if;
    if type(G, numeric) and evalb(G <= 0) then
      WARNING("null or negative shear modulus detected.");
    end if;
    if type(L, numeric) and evalb(L <= 0) then
      WARNING("null or negative length detected.");
    end if;
    if type(I_x, numeric) and evalb(I_x <= 0) then
      WARNING("null or negative x-axis inertia detected.");
    end if;
    if type(I_y, numeric) and evalb(I_y <= 0) then
      WARNING("null or negative y-axis inertia detected.");
    end if;
    if type(I_z, numeric) and evalb(I_z <= 0) then
      WARNING("null or negative z-axis inertia detected.");
    end if;
  end if;

  P_y := 12*E*I_z/(G*A*L^2);
  P_z := 12*E*I_y/(G*A*L^2);

  EA_L  := E*A/L;
  GJ_L  := G*(I_y+I_z)/L;
  EIy_L := TrussMe:-FEM:-Simplify(E*I_y/((1+P_y)*L));
  EIy_L2 := EIy_L/L; EIy_L3 := EIy_L2/L;
  EIz_L := TrussMe:-FEM:-Simplify(E*I_z/((1+P_z)*L));
  EIz_L2 := EIz_L/L; EIz_L3 := EIz_L2/L;

  out := Matrix(12, storage = sparse);
  out[1, 1]  := EA_L;       out[1, 7]   := -EA_L;
  out[2, 2]  := 12*EIz_L3;  out[2, 6]   := 6*EIz_L2;
  out[2, 8]  := -12*EIz_L3; out[2, 12]  := 6*EIz_L2;
  out[3, 3]  := 12*EIy_L3;  out[3, 5]   := -6*EIy_L2;
  out[3, 9]  := -12*EIy_L3; out[3, 11]  := -6*EIy_L2;
  out[4, 4]  := GJ_L;       out[4, 10]  := -GJ_L;
  out[5, 3]  := -6*EIy_L2;  out[5, 5]   := (4+P_y)*EIy_L;
  out[5, 9]  := 6*EIy_L2;   out[5, 11]  := (2+P_y)*EIy_L;
  out[6, 2]  := 6*EIz_L2;   out[6, 6]   := (4+P_z)*EIz_L;
  out[6, 8]  := -6*EIz_L2;  out[6, 12]  := (2+P_z)*EIz_L;
  out[7, 1]  := -EA_L;      out[7, 7]   := EA_L;
  out[8, 2]  := -12*EIz_L3; out[8, 6]   := -6*EIz_L2;
  out[8, 8]  := 12*EIz_L3;  out[8, 12]  := -6*EIz_L2;
  out[9, 3]  := -12*EIy_L3; out[9, 5]   := 6*EIy_L2;
  out[9, 9]  := 12*EIy_L3;  out[9, 11]  := 6*EIy_L2;
  out[10, 4] := -GJ_L;      out[10, 10] := GJ_L;
  out[11, 3] := -6*EIy_L2;  out[11, 5]  := (2+P_y)*EIy_L;
  out[11, 9] := 6*EIy_L2;   out[11, 11] := (4+P_y)*EIy_L;
  out[12, 2] := 6*EIz_L2;   out[12, 6]  := (2+P_z)*EIz_L;
  out[12, 8] := -6*EIz_L2;  out[12, 12] := (4+P_z)*EIz_L;
  return out;
end proc: # GetTimoshenkoBeamStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetOutputDisplacements := proc(
  fem::FEM,
  $)::Matrix;

  description "Get the output displacements of a solved and stored FEMstructure "
    "<fem>.";

  return <<fem["info"][1..-1, 1]> |
          <fem["info"][1..-1, 2]> |
          <fem["info"][1..-1, 3]> |
          fem["d"]>;
end proc: # GetOutputDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetOutputReactions := proc(
  fem::FEM,
  $)::Matrix;

  description "Get the output reactions of a solved and stored FEMstructure "
    "<fem>.";

  return <<fem["info"][1..-1, 1]> |
          <fem["info"][1..-1, 2]> |
          <fem["info"][1..-1, 3]> |
          fem["f"]>;
end proc: # GetOutputReactions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsELEMENT := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of ELEMENT type.";

  return type(var, table) and evalb(var["type"] = ELEMENT);
end proc: # IsELEMENT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsELEMENTS := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of ELEMENTS type.";

  return type(var, list(ELEMENT));
end proc: # IsELEMENTS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeElement := proc(
  name::string,
  # N1::{NODE, list({NODE, DOFS}), Vector({NODE, DOFS})},
  # N2::{NODE, list({NODE, DOFS}), Vector({NODE, DOFS})},
  # ...
  # Ni::{NODE, list({NODE, DOFS}), Vector({NODE, DOFS})},
  # K::STIFFNESS,
  {
  frame::FRAME := TrussMe:-FEM:-GenerateGenericFrame(name)
  })::ELEMENT;

  description "Make an element with name <name> on reference frame <frame>, "
    "connecting the dofs <Ni_dofs> on i-th node i <[Ni, Ni_dofs]>, connecting "
    "with stiffness <K>.";

  local nodes, dofs, i, tmp, K, L, R, Q;

  nodes := [seq(i, i = 1.._nrest-1)];
  dofs  := [seq(i, i = 1.._nrest-1)];
  for i from 1 to _nrest-1 do
    tmp := _rest[i];
    if type(tmp, NODE) then
      nodes[i] := tmp["id"];
      dofs[i]  := <0, 0, 0, 0, 0, 0>;
    elif evalb(nops(tmp) = 2) and type(tmp[1], NODE) and type(tmp[2], DOFS) then
      nodes[i] := tmp[1]["id"];
      dofs[i]  := tmp[2];
    else
      error("<N%1> must be a node or a list of node and dofs.", i);
    end if;
  end do;

  if type(_rest[-1], STIFFNESS) then
    K := _rest[-1];
  else
    error("<K> must be a stiffness matrix.");
  end if;

  return table([
    "type"      = ELEMENT,
    "name"      = name,
    "id"        = TrussMe:-FEM:-GenerateId(),
    "frame"     = frame, # Reference frame from global to local
    "nodes"     = nodes,
    "dofs"      = dofs, # Constrained dofs on nodes in local frame
    "stiffness" = K # Stiffness matrix in local frame
  ]);
end proc: # MakeElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeSpring := proc(
  name::string,
  N1::{NODE, list({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS})},
  {
  K::{algebraic, list(algebraic)} := 0,
  T::{algebraic, list(algebraic)} := 0,
  frame::FRAME                    := TrussMe:-FEM:-GenerateGenericFrame(name),
  distance::algebraic             := 0
  }, $)::ELEMENT;

  description "Make a spring element with name <name> on reference frame <frame>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with translational "
    "stiffnesses <K> or <K_x, K_y, K_z>, and torsional stiffnesses <T> or "
    "<T_x, T_y, T_z>. Optional nodes distance <distance> can be specified.";

  return TrussMe:-FEM:-MakeElement(
    name, N1, N2, TrussMe:-FEM:-GetSpringStiffness(K, T), parse("frame") = frame
  );
end proc: # MakeSpring

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeRod := proc(
  name::string,
  N1::{NODE, list({NODE, DOFS}), Vector({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS}), Vector({NODE, DOFS})},
  {
  material::MATERIAL       := TrussMe:-FEM:-MakeCarbonSteel(),
  area::algebraic          := 0,
  inertia::list(algebraic) := [0, 0, 0],
  frame::FRAME             := TrussMe:-FEM:-GenerateGenericFrame(name),
  distance::algebraic      := -1
  }, $)::ELEMENT;

  description "Make a rod element with name <name> on reference frame <frame>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with material <material>, "
    "cross-section area <area>and inertia <inertia>. Optional nodes distance "
    "<distance> can be specified.";

  local L, coordinates_1, coordinates_2;

  if evalb(distance <> -1) then
    L := distance;
  else
    if type(N1, NODE) then
      coordinates_1 := N1["coordinates"];
    elif evalb(nops(N1) = 2) and type(N1[1], NODE) and type(N1[2], DOFS) then
      coordinates_1 := N1[1]["coordinates"];
    else
      error("<N1> must be a node or a list of node and dofs.");
    end if;
    if type(N2, NODE) then
      coordinates_2 := N2["coordinates"];
    elif evalb(nops(N2) = 2) and type(N2[1], NODE) and type(N2[2], DOFS) then
      coordinates_2 := N2[1]["coordinates"];
    else
      error("<N2> must be a node or a list of node and dofs.");
    end if;
    L := TrussMe:-FEM:-Norm2(coordinates_2 - coordinates_1);
  end if;

  return TrussMe:-FEM:-MakeElement(name, N1, N2, TrussMe:-FEM:-GetRodStiffness(
      area, material["elastic_modulus"], material["shear_modulus"], L, op(inertia)
    ), parse("frame") = frame
  );
end proc: # MakeRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBeam := proc(
  name::string,
  N1::{NODE, list({NODE, DOFS}), Vector({NODE, DOFS})},
  N2::{NODE, list({NODE, DOFS}), Vector({NODE, DOFS})},
  {
  material::MATERIAL       := TrussMe:-FEM:-MakeCarbonSteel(),
  area::algebraic          := 0,
  inertia::list(algebraic) := [0, 0, 0],
  frame::FRAME             := TrussMe:-FEM:-GenerateGenericFrame(name),
  distance::algebraic      := -1,
  timoshenko::boolean      := false
  }, $)::ELEMENT;

  description "Make a beam element with name <name> on reference frame <frame>, "
    "connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting "
    "the dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with material <material> and "
    "cross-section area <area> and inertia <inertia>. Optional nodes distance "
    "<distance> and Timoshenko's beam boolean <timoshenko> can be specified.";

  local L, coordinates_1, coordinates_2;

  if evalb(distance <> -1) then
    L := distance;
  else
    if type(N1, NODE) then
      coordinates_1 := N1["coordinates"];
    elif evalb(nops(N1) = 2) and type(N1[1], NODE) and type(N1[2], DOFS) then
      coordinates_1 := N1[1]["coordinates"];
    else
      error("<N1> must be a node or a list of node and dofs.");
    end if;
    if type(N2, NODE) then
      coordinates_2 := N2["coordinates"];
    elif evalb(nops(N2) = 2) and type(N2[1], NODE) and type(N2[2], DOFS) then
      coordinates_2 := N2[1]["coordinates"];
    else
      error("<N2> must be a node or a list of node and dofs.");
    end if;
    L := TrussMe:-FEM:-Norm2(coordinates_2 - coordinates_1);
  end if;

  return TrussMe:-FEM:-MakeElement(name, N1, N2,
    `if`(timoshenko,
    TrussMe:-FEM:-GetTimoshenkoBeamStiffness(
      area, material["elastic_modulus"], material["shear_modulus"], L, op(inertia)
    ),
    TrussMe:-FEM:-GetBeamStiffness(
      area, material["elastic_modulus"], material["shear_modulus"], L, op(inertia)
    )),
    parse("frame") = frame
  );
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSTRUCTURE := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of STRUCTURE type.";

  return type(var[1], NODES) and type(var[2], ELEMENTS);
end proc: # IsSTRUCTURE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNodalDofs := proc(
  nodes::NODES,
  $)::Vector;

  description "Get nodal dofs of nodes <nodes>.";

  local dofs, i;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Getting nodal dofs... ");
  end if;

  dofs := Vector(6 * nops(nodes));
  for i from 1 to nops(nodes) do
    dofs[6*i-5..6*i] := nodes[i]["dofs"];
  end do;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  return dofs;
end proc: # GetNodalDofs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNodalDisplacements := proc(
  nodes::NODES,
  $)::Vector;

  description "Get the vector of nodal displacements of nodes <nodes>.";

  local displacements, i;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Getting nodal displacements... ");
  end if;

  displacements := Vector(6 * nops(nodes));
  for i from 1 to nops(nodes) do
    displacements[6*i-5..6*i] := nodes[i]["displacements"];
  end do;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  return displacements;
end proc: # GetNodalDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export StiffnessTransformation := proc(
  nodes::NODES,
  $)::Matrix;

  description "Compute the transformation matrix of the nodes <nodes>.";

  local T, R_transp, i;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Computing stiffness transformation matrix... ");
  end if;

  T := Matrix(6*nops(nodes), (i,j) -> `if`(evalb(i = j), 1, 0), storage = sparse);
  for i from 1 to nops(nodes) do
    R_transp := LinearAlgebra:-Transpose(TrussMe:-FEM:-Rotation(nodes[i]["frame"]));
    T[6*i-5..6*i-3, 6*i-5..6*i-3] := R_transp;
    T[6*i-2..6*i,   6*i-2..6*i  ] := R_transp;
  end do;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  return T;
end proc: # StiffnessTransformation

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GlobalStiffness := proc(
  nodes::NODES,
  elements::ELEMENTS,
  $)::Matrix;

  description "Compute the global stiffness matrix of the nodes <nodes> and "
    "elements <elements>.";

  local K, Q_i, D_i, K_i, R_i, n_i, p_i, d_i, i, j, k;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Computing global stiffness matrix... ");
  end if;

  K := Matrix(6*nops(nodes), storage = sparse);
  for i from 1 to nops(elements) do

    # Retrieve number of nodes
    n_i := nops(elements[i]["nodes"]);

    # Build rotation matrix Q, retrieve nodes positions and build diagonal
    Q_i := Matrix(6*n_i, storage = sparse);
    R_i := TrussMe:-FEM:-Rotation(elements[i]["frame"]);
    p_i := Vector(n_i);
    d_i := Vector(6*n_i);
    for j from 1 to n_i do
      Q_i[6*j-5..6*j-3, 6*j-5..6*j-3] := R_i;
      Q_i[6*j-2..6*j,   6*j-2..6*j  ] := R_i;
      p_i[j] := TrussMe:-FEM:-GetObjById(
        nodes, elements[i]["nodes"][j], parse("position") = true
      );
      d_i[6*j-5..6*j] := 1 -~ elements[i]["dofs"][j];
    end do;

    # Element stiffness contribution selecting only constrained dofs (= 0)
    D_i := LinearAlgebra:-DiagonalMatrix(d_i);
    K_i := Q_i.(D_i.elements[i]["stiffness"].D_i).LinearAlgebra:-Transpose(Q_i);

    # Insert element stiffness contribution in global stiffness matrix
    for j from 1 to n_i do
      for k from 1 to n_i do
        K[6*p_i[j]-5..6*p_i[j], 6*p_i[k]-5..6*p_i[k]] :=
        K[6*p_i[j]-5..6*p_i[j], 6*p_i[k]-5..6*p_i[k]] + K_i[6*j-5..6*j, 6*k-5..6*k];
      end do;
    end do;
  end do;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  return K;
end proc: # GlobalStiffness

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GlobalStiffnessPrime := proc(
  nodes::NODES,
  elements::ELEMENTS,
  $)::Matrix;

  description "Compute the global stiffness matrix of the nodes <nodes> and "
    "elements <elements>.";

  local T, K;

  T := TrussMe:-FEM:-StiffnessTransformation(nodes);
  K := TrussMe:-FEM:-GlobalStiffness(nodes, elements);
  return T.K.LinearAlgebra:-Transpose(T);
end proc: # GlobalStiffnessPrime

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsFEM := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of FEM type.";

  return type(var, table) and evalb(var["type"] = FEM);
end proc: # IsFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SplitFEM := proc(
  fem::FEM,
  $)

  description "Split the FEM structure <fem> in free and constrained dofs.";

  local d, K, f, n, i;

  # Get permutation and unpermutation vectors
  fem["perm"]   := sort(fem["dofs"], `>`, output = 'permutation');
  fem["unperm"] := [seq(i, i = 1..nops(fem["perm"]))];
  for i from 1 to nops(fem["unperm"]) do
    fem["unperm"][fem["perm"][i]] := i;
  end do;

  # Compute global stiffness matrix, displacements and loads vectors
  fem["K"] := TrussMe:-FEM:-Simplify(
    TrussMe:-FEM:-GlobalStiffnessPrime(fem["nodes"], fem["elements"])
  );
  fem["d"] := TrussMe:-FEM:-Simplify(
    TrussMe:-FEM:-GetNodalDisplacements(fem["nodes"])
  );
  fem["f"] := TrussMe:-FEM:-Simplify(
    TrussMe:-FEM:-GetNodalLoads(fem["nodes"], fem["loads"])
  );

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Splitting system... ");
  end if;

  # Compute permuted system stiffness matrices, displacements and loads vectors
  K := fem["K"][fem["perm"], fem["perm"]];
  d := fem["d"][fem["perm"]];
  f := fem["f"][fem["perm"]];
  n := add(fem["dofs"]);
  fem["info_f"] := fem["info"][fem["perm"]][1..n, 1..-1];
  fem["info_s"] := fem["info"][fem["perm"]][n+1..-1, 1..-1];
  fem["K_ff"]   := K[1..n, 1..n];
  fem["K_fs"]   := K[1..n, n+1..-1];
  fem["K_sf"]   := K[n+1..-1, 1..n];
  fem["K_ss"]   := K[n+1..-1, n+1..-1];
  fem["d_f"]    := Vector(0);
  fem["d_s"]    := d[n+1..-1];
  fem["f_f"]    := f[1..n];
  fem["f_s"]    := Vector(0);
  fem["f_r"]    := f[n+1..-1];

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  # Veiling variables
  fem["label"]  := "";
  fem["veils"]  := [];
  fem["solved"] := false;
  return NULL;
end proc: # SplitFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export RecastFEM := proc(
  fem::FEM,
  $)

  description "Recast the system <fem> to avoid singularities.";

  local i, fill_tmp, info_tmp, idx_tmp;

  # Get system fill-in matrix and info
  fill_tmp := map(x -> `if`(x = 0, 0, 1), fem["K_ff"]);
  info_tmp := fem["info"][fem["perm"]];
  idx_tmp  := [];

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Looking for singularities... ");
  end if;

  # Check for empty rows
  for i from 1 to LinearAlgebra:-RowDimension(fem["K_ff"]) do
    if evalb(add(fill_tmp[i, 1..-1]) = 0) then
      if evalb(fem["f_f"][i] = 0) then
        idx_tmp := [op(idx_tmp), i];
        if TrussMe:-FEM:-m_WarningMode then
          WARNING("recasting system at node, dof = [name: %1, id: %2, direction: "
            "%3] is unconstrained and has a zero load.", op(info_tmp[i][1..3]));
        end if;
      else
        error("singular system detected, dof = [name: %1, id: %2, direction: %3] "
          "is unconstrained and has a non-zero load.", op(info_tmp[i][1..3]));
      end if;
    end if;
  end do;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  # Recast nodal dofs and split fem structure again
  for i in idx_tmp do
    fem["dofs"][fem["perm"][i]] := 0;
  end do;
  TrussMe:-FEM:-SplitFEM(fem);
  return NULL;
end proc: # RecastFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateFEM := proc(
  nodes::NODES,
  elements::ELEMENTS,
  loads::LOADS,
  {
  tryhard::boolean := false
  }, $)::FEM;

  description "Generate a FEM structure from the nodes <nodes>, elements "
    "<elements>, and loads <loads>. If <tryhard> is enabled, the stiffness "
    "matrix will be recasted to avoid singularities.";

  local fem, direction;

  # Create FEM structure
  fem := table([
    "type"     = FEM,
    "nodes"    = nodes,
    "elements" = elements,
    "loads"    = loads,
    "info"     = map(
      x -> seq([x["name"], x["id"], direction],
      direction = ["tx", "ty", "tz", "rx", "ry", "rz"]
    ), nodes)
  ]);

  # Split free and constrained dofs
  fem["dofs"] := TrussMe:-FEM:-GetNodalDofs(fem["nodes"]);
  TrussMe:-FEM:-SplitFEM(fem);

  # Recast the stiffness matrix to avoid singularities
  if tryhard then
    TrussMe:-FEM:-RecastFEM(fem);
  end if;

  # Return the FEM structure
  return fem;
end proc: # GenerateFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SolveFEM := proc(
  fem::FEM,
  {
  use_LAST::boolean     := false,
  use_LEM::boolean      := true,
  use_SIG::boolean      := true,
  factorization::string := "LU",
  label::string         := "V"
  }, $)

  description "Solve the FEM structure <fem> and optionally use LAST LU "
    "decompostion <use_LAST> and veil the expressions <use_LEM> with signature "
    "mode <use_SIG> and label <label>. Factorization method <factorization> "
    "can be choosen between 'LU', fraction-free 'FFLU', 'QR', and Gauss-Jordan "
    "'GJ'.";

  local dofs, LAST_obj, LEM_obj, i;

  # Solve free dofs
  if use_LAST then
    try
      if TrussMe:-FEM:-m_VerboseMode then
        printf("Initializing LAST and LEM objects... ");
      end if;
      LAST_obj := Object(LAST);
      LAST_obj:-InitLEM(LAST_obj, label);
      LEM_obj := LAST_obj:-GetLEM(LAST_obj);
      if TrussMe:-FEM:-m_VerboseMode then
        printf("DONE\n");
      end if;
    catch:
      if TrussMe:-FEM:-m_VerboseMode then
        printf(" FAILED\n");
      end if;
      error("LAST or LEM package not installed, please refere to the "
        "documentation for more information.");
    end try;

    # Set verbose and warning modes
    LAST_obj:-SetVerboseMode(LAST_obj, TrussMe:-FEM:-m_VerboseMode);
    LAST_obj:-SetWarningMode(LAST_obj, TrussMe:-FEM:-m_WarningMode);

    # Set signature checking and time limit
    LEM_obj:-SetSignatureMode(LEM_obj, use_SIG);
    LAST_obj:-SetTimeLimit(LAST_obj, TrussMe:-FEM:-m_TimeLimit);

    if TrussMe:-FEM:-m_VerboseMode then
      printf("Performing matrix factorization... ");
    end if;

    # Perform decomposition
    if evalb(factorization = "LU") then
      LAST_obj:-LU(LAST_obj, fem["K_ff"]);
    elif evalb(factorization = "FFLU") then
      LAST_obj:-FFLU(LAST_obj, fem["K_ff"]);
    elif evalb(factorization = "QR") then
      LAST_obj:-QR(LAST_obj, fem["K_ff"]);
    elif evalb(factorization = "GJ") then
      LAST_obj:-GJ(LAST_obj, fem["K_ff"]);
    else
      error("<factorization> must be 'LU', 'FFLU', 'QR', or 'GJ'.");
    end if;

    if TrussMe:-FEM:-m_VerboseMode then
      printf("DONE\n");
      printf("Solving system deformations... ");
    end if;

    # Solve deformations
    fem["d_f"] := LAST_obj:-SolveLinearSystem(
      LAST_obj, fem["f_f"] - fem["K_fs"].fem["d_s"]
    );

    # Unveil expressions if required
    if use_LEM then
      fem["label"] := label;
      fem["veils"] := LEM_obj:-VeilList(LEM_obj);
    else
      fem["d_f"]   := LEM_obj:-Unveil(LEM_obj, fem["d_f"]);
      fem["label"] := "";
      fem["veils"] := [];
    end if;
  else

    if TrussMe:-FEM:-m_VerboseMode then
      printf("Solving system deformations... ");
    end if;

    # Solve deformations
    fem["d_f"] := LinearAlgebra:-LinearSolve(
      fem["K_ff"], fem["f_f"] - fem["K_fs"].fem["d_s"]
    ) assuming real;
    fem["label"] := "";
    fem["veils"] := [];
  end if;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
    printf("Solving system reactions... ");
  end if;

  # Solve reactions
  fem["f_s"] := fem["K_sf"].fem["d_f"] + fem["K_ss"].fem["d_s"] - fem["f_r"];

  # Store output displacements and reactions and restore initial permutation
  fem["d"] := convert(<fem["d_f"], fem["d_s"]>[fem["unperm"]], Vector);
  fem["f"] := convert(<fem["f_f"], fem["f_s"]>[fem["unperm"]], Vector);

  # Mark as solved
  fem["solved"] := true;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  return NULL;
end proc: # SolveFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export StoreFEM := proc(
  fem::FEM,
  $)

  description "Store the FEM structure <fem> in the nodes <nodes>.";

  local i;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Storing system... ");
  end if;

  for i from 1 to nops(fem["nodes"]) do
    fem["nodes"][i]["output_displacements"] := fem["d"][6*i-5..6*i];
    fem["nodes"][i]["output_reactions"]     := fem["f"][6*i-5..6*i];
  end do;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;
  return NULL;
end proc: # StoreFEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
