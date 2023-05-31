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

  return type(var, LOAD);
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
  RF::FRAME,
  node::NODE,
  components::COMPONENTS,
  $)::LOAD;

  description "Create a load with name <name> acting on the node (or its id) "
    "<node> on the frame <RF> with components <components>.";

  if evalb(nops(components) <> 6) then
    error("<displacements> must be a list of 6 elements.");
  elif evalb(add(-(dofs-2) *~ components) <> 0) then
    error("<displacements> must be defined only for constrained dofs.");
  end if;

  return table(
    "type"       = LOAD,
    "name"       = name,
    "id"         = TrussMe_FEM:-GenerateId(),
    "frame"      = RF,
    "node"       = node["id"],
    "components" = components
  );
end proc: # MakeLoad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsLOADS::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of LOADS type.";

  return type(var, list(LOAD));
end proc: # IsLOADS

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetLocalLoads::static := proc(
  RF::FRAME,
  node::NODE,
  $)::LOADS;

  description "Get the loads acting on the node (or its id) <node> on the "
    "frame <RF>.";

end proc: # GetLocalLoads

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
