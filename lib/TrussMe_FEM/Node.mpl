# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           _   _           _                                 #
#                          | \ | | ___   __| | ___                            #
#                          |  \| |/ _ \ / _` |/ _ \                           #
#                          | |\  | (_) | (_| |  __/                           #
#                          |_| \_|\___/ \__,_|\___|                           #
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

  return type(var, Vector) and evalb(LinearAlgebra:-Dimension(var) = 6) and
    not has(map(has, var, {0, 1}), false);
end proc: # IsDOFS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeNode::static := proc(
  name::string,
  coords::POINT,
  RF::FRAME := Matrix(4, shape = identity),
  dofs::DOFS := <1, 1, 1, 1, 1, 1>,
  $)::NODE;

  description "Create a node with name <name> and coordinates <coords> in the "
    "reference frame <RF>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained to the ground "
    "in the direction given from <RF>."

  if evalb(LinearAlgebra:-Dimension(dofs) = 6) or
     has(map(has, dofs, {0, 1}), false) then
    error("dofs must have 6 elements, each of which is either 0 or 1.");
  end if;

  return table([
    "type"          = NODE,
    "coordinates"   = coords,
    "name"          = name,
    "id"            = TrussMe_FEM:-GenerateId(),
    "frame"         = RF,
    "dofs"          = dofs,
    "loads"         = [],
    "displacements" = []
  ]);
end proc: # MakeNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCompliantNode := proc(
  name::string,
  coords::{list, POINT, VECTOR},
  RF::FRAME := Matrix(4, shape = identity),
  dofs::DOFS := <1, 1, 1, 1, 1, 1>,
  {
  K::{algebraic, list(algebraic)},
  T::{algebraic, list(algebraic)}
  }, $)::NODE, ELEMENT, NODE;

  description "Create a node with name <name> and coordinates <coords> in the "
    "reference frame <RF>. The constraints on dofs are specified by <dofs>, "
    "where 1 means free and 0 means that the dof is constrained to the ground "
    "in the direction given from <RF>. The node is also connected to a "
    "compliant spring element with traslational stiffness <K> and torsional "
    "stiffness <T>.";

  local L, R, C;

  L := TrussMe_FEM:-MakeNode(name, coords, RF, dofs);
  R := TrussMe_FEM:-MakeNode(cat(name, "_compliant"), coords, RF);
  C := TrussMe_FEM:-MakeSpring(
    cat(name, "_spring"), parse("K") = K, parse("T") = T
  );
  return L, C, R;
end proc: # MakeCompliantNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
