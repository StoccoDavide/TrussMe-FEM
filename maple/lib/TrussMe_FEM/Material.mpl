# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                     __  __       _            _       _                     #
#                    |  \/  | __ _| |_ ___ _ __(_) __ _| |                    #
#                    | |\/| |/ _` | __/ _ \ '__| |/ _` | |                    #
#                    | |  | | (_| | ||  __/ |  | | (_| | |                    #
#                    |_|  |_|\__,_|\__\___|_|  |_|\__,_|_|                    #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsMATERIAL := proc(
  var::anything,
  $)::boolean;

  description "Check if the variable <var> is of MATERIAL type.";

  return type(var, table) and evalb(var["type"] = MATERIAL);
end proc: # IsMATERIAL

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeMaterial := proc({
  name::string               := "Undefined",
  elastic_modulus::algebraic := 0,
  poisson_ratio::algebraic   := 0,
  shear_modulus::algebraic   := elastic_modulus/(2*(1+poisson_ratio)),
  density::algebraic         := 0
  }, $)::MATERIAL;

  description "Define material with inputs: name <name>, elastic modulus (E) "
    "<elastic_modulus>, Poisson's ratio (nu) <poisson_ratio>, shear modulus "
    "<shear_modulus> (default = E/(2*(1+nu))), and density <density>.";

  return table([
    "type"            = MATERIAL,
    "name"            = name,
    "id"              = TrussMe_FEM:-GenerateId(),
    "elastic_modulus" = elastic_modulus,
    "shear_modulus"   = shear_modulus,
    "poisson_ratio"   = poisson_ratio,
    "density"         = density
  ]);
end proc: # MakeMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCarbonSteel := proc( $ )::MATERIAL;

  description "Get default steel material with 'CarbonSteel' name, elastic "
    "modulus E = 210.0e+09 (Pa), Poisson's ratio nu = 0.3 (-), shear modulus"
    "E/(2*(1+nu), density rho = 7850.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "CarbonSteel",
    parse("elastic_modulus") = 210.0e+09,
    parse("poisson_ratio")   = 0.3,
    parse("density")         = 7850.0
  );
end proc: # CarbonSteel

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeInoxSteel := proc( $ )::MATERIAL;

    description "Get default steel material with 'InoxSteel' name, elastic "
      "modulus E = 200.0e+09 (Pa), Poisson's ratio nu = 0.3 (-), shear modulus"
      "E/(2*(1+nu), density rho = 8000.0 (kg/m^3).";

    return TrussMe_FEM:-MakeMaterial(
      parse("name")            = "InoxSteel",
      parse("elastic_modulus") = 200.0e+09,
      parse("poisson_ratio")   = 0.3,
      parse("density")         = 8000.0
    );
  end proc: # InoxSteel

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeTitanium := proc( $ )::MATERIAL;

  description "Get default titanium material with 'Titanium' name, elastic "
    "modulus E = 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 4500.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Titanium",
    parse("elastic_modulus") = 110.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 4500.0
  );
end proc: # Titanium

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCopper := proc( $ )::MATERIAL;

  description "Get default copper material with 'Copper' name, elastic "
    "modulus E = 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 4500.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Copper",
    parse("elastic_modulus") = 110.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 4500.0
  );
end proc: # Copper

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBrass := proc( $ )::MATERIAL;

  description "Get default brass material with 'Brass' name, elastic "
    "modulus E = 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 4500.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Brass",
    parse("elastic_modulus") = 110.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 4500.0
  );
end proc: # Brass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBronze := proc( $ )::MATERIAL;

  description "Get default bronze material with 'Bronze' name, elastic "
    "modulus E = 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 4500.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Bronze",
    parse("elastic_modulus") = 110.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 4500.0
  );
end proc: # Bronze

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeLead := proc( $ )::MATERIAL;

  description "Get default lead material with 'Lead' name, elastic "
    "modulus E = 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 4500.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Lead",
    parse("elastic_modulus") = 110.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 4500.0
  );
end proc: # Lead

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeZinc := proc( $ )::MATERIAL;

  description "Get default zinc material with 'Zinc' name, elastic "
    "modulus E = 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 4500.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Zinc",
    parse("elastic_modulus") = 110.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 4500.0
  );
end proc: # Zinc

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeMagnesium  := proc( $ )::MATERIAL;

  description "Get default magnesium material with 'Magnesium' name, elastic "
    "modulus E = 45.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 1800.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Magnesium",
    parse("elastic_modulus") = 45.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 1800.0
  );
end proc: # Magnesium

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeAlluminium := proc( $ )::MATERIAL;

  description "Get default alluminium material with 'Alluminium' name, elastic "
    "modulus E = 69.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 8000.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Alluminium",
    parse("elastic_modulus") = 70.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 8000.0
  );
end proc: # Alluminium

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeAvional := proc( $ )::MATERIAL;

  description "Get default avional material with 'Avional' name, elastic "
    "modulus E = 70.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 2690.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Avional",
    parse("elastic_modulus") = 70.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 2690.0
  );
end proc: # Avional

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakePeraluman := proc( $ )::MATERIAL;

  description "Get default peraluman material with 'Peraluman' name, elastic "
    "modulus E = 70.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 2700.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Peraluman",
    parse("elastic_modulus") = 70.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 2700.0
  );
end proc: # Peraluman

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeAnticorodal := proc( $ )::MATERIAL;

  description "Get default anticorodal material with 'Anticorodal' name, elastic "
    "modulus E = 69.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 2700.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Anticorodal",
    parse("elastic_modulus") = 69.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 2700.0
  );
end proc: # Anticorodal

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeCarpental := proc( $ )::MATERIAL;

  description "Get default carpental material with 'Carpental' name, elastic "
    "modulus E = 72.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 2780.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Carpental",
    parse("elastic_modulus") = 72.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 2780.0
  );
end proc: # Carpental

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeErgal := proc( $ )::MATERIAL;

  description "Get default ergal material with 'Ergal' name, elastic "
    "modulus E = 72.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulus"
    "E/(2*(1+nu), density rho = 2780.0 (kg/m^3).";

  return TrussMe_FEM:-MakeMaterial(
    parse("name")            = "Ergal",
    parse("elastic_modulus") = 72.0e+09,
    parse("poisson_ratio")   = 0.33,
    parse("density")         = 2780.0
  );
end proc: # Ergal

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
