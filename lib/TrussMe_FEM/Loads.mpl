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

export IsFORCE::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of FORCE type.";

  return type(var, table) and evalb(var["type"] = FORCE);
end proc: # IsFORCE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsQFORCE::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of QFORCE type.";

  return type(var, table) and evalb(var["type"] = QFORCE);
end proc: # IsQFORCE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsMOMENT::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of MOMENT type.";

  return type(var, table) and evalb(var["type"] = MOMENT);
end proc: # IsMOMENT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsQMOMENT::static := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of QMOMENT type.";

  return type(var, table) and evalb(var["type"] = QMOMENT);
end proc: # IsQMOMENT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
