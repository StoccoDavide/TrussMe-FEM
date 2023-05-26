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

export MakeNode::static := proc(
  name::string,
  coords::POINT,
  RF::FRAME,
  $)::NODE;

  return table([
    "type"          = NODE,
    "coordinates"   = coords,
    "name"          = name,
    "id"            = TrussMe_FEM:-GenerateId(),
    "frame"         = RF,
    "loads"         = [],
    "displacements" = [],
    "reactions"     = []
  ]);
end proc: # MakeNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCompliantNode := proc(
  name::string,
  coords::{list, POINT, VECTOR},
  constrained_dof::list,
  RF::FRAME,
  K::{algebraic, list(algebraic)},
  T::{algebraic, list(algebraic)}
  $)::NODE, ELEMENT, NODE;

  description

  local L, R, C;

  L := TrussMe_FEM:-MakeNode(name, coords, RF);
  R := TrussMe_FEM:-MakeNode(cat(name, "_compliant"), coords, RF);
  C := TrussMe_FEM:-MakeSpring(cat(name, "_spring"), K, T)
  return L, C, R;
end proc: # MakeCompliantNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

