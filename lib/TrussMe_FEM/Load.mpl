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

export IsLOAD := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of LOAD type.";

  return type(var, table) and var["type"] = LOAD;
end proc: # IsLOAD

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsCOMPONENTS := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of COMPONENTS type.";

  return type(var, Vector) and evalb(LinearAlgebra:-Dimension(var) = 6);
end proc: # IsCOMPONENTS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeLoad := proc(
  name::string,
  node::NODE,
  components::COMPONENTS,
  {
  frame::{string, FRAME} := node["id"]
  }, $)::LOAD;

  description "Create a load with name <name> acting on the node (or its id) "
    "<node> on the frame <frame> with components <components>.";
print("aa");
  if evalb(add(-(node["dofs"] -~ 1) *~ components) <> 0) then
    error("components must be defined only for constrained dofs.");
  end if;
print("aa");

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

export IsLOADS := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of LOADS type.";

  return type(var, list(LOAD));
end proc: # IsLOADS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNodalLoads := proc(
  nodes::NODES,
  loads::LOADS,
  $)::Vector;

  description "Get the vector of nodal loads <loads> of nodes <nodes>.";

  local F, i, j, R;

  F := Vector(6 * nops(nodes), storage = sparse);
  for i from 1 to nops(loads) do
    # Node position
    j := TrussMe_FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
    # Node loads in node frame
    if evalb(loads[i]["frame"] = nodes[j]["id"]) then
      R := Matrix(3, shape = identity);
    else
      R := (LinearAlgebra:-Transpose(nodes[j]["frame"]).loads[i]["frame"])[1..3, 1..3];
    end if;
    F[6*j-5..6*j] := <
      R.loads[i]["components"][1..3], R.loads[i]["components"][4..6]
    >;
  end do;
  return F;
end proc: # GetNodalLoads

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
