# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                         _                    _                              #
#                        | |    ___   __ _  __| |___                          #
#                        | |   / _ \ / _` |/ _` / __|                         #
#                        | |__| (_) | (_| | (_| \__ \                         #
#                        |_____\___/ \__,_|\__,_|___/                         #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsLOAD::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of LOAD type.";

  return type(var, table) and var["type"] = LOAD;
end proc: # IsLOAD

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsCOMPONENTS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of COMPONENTS type.";

  return type(var, list({algebraic, procedure})) and evalb(nops(var) = 6);
end proc: # IsCOMPONENTS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeLoad::static := proc(
  name::string,
  node::NODE,
  components::COMPONENTS,
  {
  frame::{string, FRAME} := node["id"]
  }, $)::LOAD;

  description "Create a load with name <name> acting on the node (or its id) "
    "<node> on the frame <frame> with components <components>.";

  if evalb(nops(components) <> 6) then
    error("<displacements> must be a list of 6 elements.");
  elif evalb(add(-(node["dofs"] -~ 1) *~ components) <> 0) then
    error("<displacements> must be defined only for constrained dofs.");
  end if;

  return table([
    "type"       = LOAD,
    "name"       = name,
    "id"         = TrussMe_FEM:-GenerateId(),
    "frame"      = frame,
    "node"       = node["id"],
    "components" = components
  ]);
end proc: # MakeLoad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsLOADS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of LOADS type.";

  return type(var, list(LOAD));
end proc: # IsLOADS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNodalLoads::static := proc(
  nodes::NODES,
  loads::LOADS,
  $)::Vector;

  description "Get the vector of nodal loads <loads> of nodes <nodes>.";

  local F, i, j, R;

  printf("TrussMe_FEM:-GetNodalLoads() ... ");
  F := Vector(6 * nops(nodes), storage = sparse);
  for i from 1 to nops(loads) do
    # Node position
    j := TrussMe_FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
    # Node loads in node frame
    if evalb(loads[i]["frame"] = nodes[j]["id"]) then
      R := LinearAlgebra:-IdentityMatrix(3);
    else
      R := (LinearAlgebra:-Transpose(nodes[j]["frame"]).loads[i]["frame"])[1..3, 1..3];
    end if;
    F[6*j-5..6*j] := <
      R.<op(loads[i]["components"][1..3])>, R.<op(loads[i]["components"][4..6])>
    >;
  end do;
  printf("DONE\n");
  return convert(F, Vector);
end proc: # GetNodalLoads

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
