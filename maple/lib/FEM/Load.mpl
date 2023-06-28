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

  return table([
    "type"       = LOAD,
    "name"       = name,
    "id"         = TrussMe:-FEM:-GenerateId(),
    "frame"      = frame, # Reference frame from global to local or node id
    "node"       = node["id"],
    "components" = components # Components in the local frame
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

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Getting nodal loads... ");
  end if;

  F := Vector(6 * nops(nodes), storage = sparse);
  for i from 1 to nops(loads) do
    # Node position
    j := TrussMe:-FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
    # Node loads in node frame
    if evalb(loads[i]["frame"] = nodes[j]["id"]) then
      R := Matrix(3, shape = identity);
    else
      R := TrussMe:-FEM:-Rotation(TrussMe:-FEM:-InverseFrame(nodes[j]["frame"]).loads[i]["frame"]);
    end if;
    F[6*j-5..6*j] := F[6*j-5..6*j] + convert(
      <R.loads[i]["components"][1..3], R.loads[i]["components"][4..6]>, Vector
    );
  end do;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  return F;
end proc: # GetNodalLoads

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
