# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#           _____                   __  __        _____ _____ __  __          #
#          |_   _| __ _   _ ___ ___|  \/  | ___  |  ___| ____|  \/  |         #
#            | || '__| | | / __/ __| |\/| |/ _ \ | |_  |  _| | |\/| |         #
#            | || |  | |_| \__ \__ \ |  | |  __/ |  _| | |___| |  | |         #
#            |_||_|   \__,_|___/___/_|  |_|\___| |_|   |_____|_|  |_|         #
#              A Maple package for symbolic FEM structural analysis           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

unprotect('TrussMe_FEM');
TrussMe_FEM := module()

  description "A Maple package for symbolic FEM structural analysis.";

  option object;

  local m_ground                := NULL;
  local m_gravity               := NULL;
  local m_earth                 := NULL;
  local m_VerboseMode           := false;
  local m_WarningMode           := true;
  local m_TimeLimit             := 5;
  local m_StoredData            := [];
  local m_BeamColor             := "SteelBlue";
  local m_RodColor              := "Niagara DarkOrchid";
  local m_RigidBodyColor        := "Indigo";
  local m_CompliantSupportColor := "DarkGreen";
  local m_SupportColor          := "DarkOrange";
  local m_CompliantJointColor   := "LightSalmon";
  local m_JointColor            := "MediumSeaGreen";
  local m_EarthColor            := "Firebrick";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print 'TrussMe_FEM' module information.";

    printf(
      "+------------------------------------------------------------------------------+\n"
      "| 'TrussMe_FEM' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023 |\n"
      "| Current version authors:                                                     |\n"
      "|   Matteo Larcher and Davide Stocco.                                          |\n"
      "+------------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad::static := proc()

    description "'TrussMe_FEM' module load procedure.";

    local lib_base_path, i, types;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("TrussMe_FEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error "cannot find 'TrussMe_FEM' library.";
    end if;

    # Define module types
    TypeTools:-AddType('FRAME',  TrussMe_FEM:-IsFrame);  protect('FRAME');
    TypeTools:-AddType('VECTOR', TrussMe_FEM:-IsVector); protect('VECTOR');
    TypeTools:-AddType('EARTH',  TrussMe_FEM:-IsEarth);  protect('EARTH');
    TypeTools:-AddType('NODE',   TrussMe_FEM:-IsNode);   protect('NODE');
    #TypeTools:-AddType('ELEMENT', TrussMe_FEM:-IsElement); protect('ELEMENT')
    #TypeTools:-AddType('FORCE', TrussMe_FEM:-IsForce); protect('FORCE')
    #TypeTools:-AddType('MOMENT', TrussMe_FEM:-IsMoment); protect('MOMENT')
    #TypeTools:-AddType('QFORCE', TrussMe_FEM:-IsQForce); protect('QFORCE')
    #TypeTools:-AddType('QMOMENT', TrussMe_FEM:-IsQMoment); protect('QMOMENT')
    #TypeTools:-AddType('SUPPORT', TrussMe_FEM:-IsSupport); protect('SUPPORT')
    #TypeTools:-AddType('JOINT',   TrussMe_FEM:-IsJoint); protect('JOINT')
    #TypeTools:-AddType('MATERIAL', TrussMe_FEM:-IsMaterial); protect('MATERIAL')
    #TypeTools:-AddType('STRUCTURE', TrussMe_FEM:-IsStructure); protect('STRUCTURE')
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()

    description "TrussMe_FEM' module unload procedure.";

    unprotect(
      'EARTH',
      'FRAME',
      'POINT',
      'VECTOR',
      'BEAM',
      'ROD',
      'RIGID_BODY',
      'FORCE',
      'MOMENT',
      'QFORCE',
      'QMOMENT',
      'SUPPORT',
      'JOINT',
      'MATERIAL',
      'STRUCTURE'
    );
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleCopy::static := proc(
    _self::TrussMe_FEM,
    proto::TrussMe_FEM,
    $)

    description "Copy the 'TrussMe_FEM' object <proto> into <self>.";


  end proc: # ModuleCopy


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Protect::static := proc()

    description "Protect module type symbols.";

    protect(
      'EARTH',
      'FRAME',
      'POINT',
      'VECTOR',
      'BEAM',
      'ROD',
      'RIGID_BODY',
      'FORCE',
      'MOMENT',
      'QFORCE',
      'QMOMENT',
      'SUPPORT',
      'JOINT',
      'MATERIAL',
      'STRUCTURE'
    );
    return NULL;
  end proc: # Protect

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Unprotect::static := proc()

    description "Unprotect module type symbols.";

    unprotect(
      'EARTH',
      'FRAME',
      'POINT',
      'VECTOR',
      'BEAM',
      'ROD',
      'RIGID_BODY',
      'FORCE',
      'MOMENT',
      'QFORCE',
      'QMOMENT',
      'SUPPORT',
      'JOINT',
      'MATERIAL',
      'STRUCTURE'
    );
    return NULL;
  end proc: # Unprotect

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Initialize::static := proc(
    _self::TrussMe_FEM,
    {
    gruond::FRAME   := Matrix(4, shape = identity),
    gravity::VECTOR := Vector(3)
    }, $)

    description "TrussMe_FEM' module initialization.";

    _self:-m_ground  := ground;
    _self:-m_gravity := gravity;
    _self:-m_earth   := table([
      "type"   = EARTH,
      "name"   = "earth",
      "frame"  = ground
    ]):
    return NULL;
  end proc: # Initialize

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetOptions::static := proc(
    _self::TrussMe_FEM,
    {
      VerboseMode::{boolean, nothing}   := NULL,
      WarningMode::{boolean, nothing}   := NULL,
      TimeLimit::{nonnegative, nothing} := NULL
    }, $)

    description "Set the module options: verbose mode <VerboseMode>, warning "
      "mode <WarningMode>, and time limit <TimeLimit>.";

    if (VerboseMode <> NULL) then
      _self:-SetVerboseMode(_self, VerboseMode);
    end if;

    if (WarningMode <> NULL) then
      _self:-SetWarningMode(_self, WarningMode);
    end if;

    if (TimeLimit <> NULL) then
      _self:-SetTimeLimit(_self, TimeLimit);
    end if;
    return NULL;
  end proc: # SetModuleOptions

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerboseMode::static := proc(
    _self::TrussMe_FEM,
    $)

    description "Enable the verbose mode of the module.";

    _self:-m_VerboseMode := true;
    return NULL;
  end proc: # EnableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerboseMode::static := proc(
    _self::TrussMe_FEM,
    $)

    description "Disable the verbose mode of the module.";

    _self:-m_VerboseMode := false;
    return NULL;
  end proc: # DisableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetVerboseMode::static := proc(
    _self::TrussMe_FEM,
    $)::boolean;

    description "Get the verbose mode of the module.";

    return _self:-m_VerboseMode;
  end proc: # GetVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVerboseMode::static := proc(
    _self::TrussMe_FEM,
    verbose_mode::boolean,
    $)

    description "Set the verbose mode of the module to <verbose_mode>.";

    _self:-m_VerboseMode := verbose_mode;
    return NULL;
  end proc: # SetVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableWarningMode::static := proc(
    _self::TrussMe_FEM,
    $)

    description "Enable the warning mode of the module.";

    _self:-m_WarningMode := true;
    return NULL;
  end proc: # EnableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableWarningMode::static := proc(
    _self::TrussMe_FEM,
    $)

    description "Disable the warning mode of the module.";

    _self:-m_WarningMode := false;
    return NULL;
  end proc: # DisableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetWarningMode::static := proc(
    _self::TrussMe_FEM,
    $)::boolean;

    description "Get the warning mode of the module.";

    return _self:-m_WarningMode;
  end proc: # GetWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWarningMode::static := proc(
    _self::TrussMe_FEM,
    warning_mode::boolean,
    $)

    description "Set the warning mode of the module to <warning_mode>.";

    _self:-m_WarningMode := warning_mode;
    return NULL;
  end proc: # SetWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetTimeLimit::static := proc(
    _self::TrussMe_FEM,
    $)::nonnegative;

    description "Get the time limit of the module.";

    return _self:-m_TimeLimit;
  end proc: # GetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetTimeLimit::static := proc(
    _self::TrussMe_FEM,
    time_limit::nonnegative,
    $)

    description "Set the time limit of the module to <time_limit>.";

    _self:-m_TimeLimit := time_limit;
    return NULL;
  end proc: # SetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetGround::static := proc(
    _self::TrussMe_FEM,
    ground::FRAME,
    $)

    description "Set ground reference frame to <ground>.";

    _self:-m_ground := ground;
    if (_self:-WarningMode) then
      WARNING("ground reference is changed.");
    end if;
    _self:-m_earth["frame"] := ground;
    return NULL;
  end proc: # SetGround

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetGround::static := proc(
    _self::TrussMe_FEM,
    $)::FRAME;

    description "Get ground reference frame.";

    return _self:-m_ground;
  end proc: # GetGround

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetGravity::static := proc(
    _self::TrussMe_FEM,
    g::VECTOR,
    $)

    description "Set gravity vector with [x, y, z]^T components of <vec>.";

    if evalb(LinearAlgebra:-Dimension(g) = 3) then
      m_gravity := g;
    else
      error "invalid gravity vector dimension.";
    end if;
    return NULL;
  end proc: # SetGravity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetGravity::static := proc(
    _self::TrussMe_FEM,
    $)::VECTOR;

    description "Get gravity vector.";

    return _self:-m_gravity;
  end proc: # GetGravity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsEarth::static := proc(
    var::anything,
    $)::boolean;

    description "Check if the variable <var> is of EARTH type.";

    return type(var, table) and evalb(var["type"] = EARTH);
  end proc: # IsEarth

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeMaterial::static := proc({
      name::string               := "DeafultSteel",
      elastic_modulus::algebraic := 210.0e+09,
      poisson_ratio::algebraic   := 0.3,
      shear_modulus::algebraic   := elastic_modulus/(2*(1+poisson_ratio)),
      density::algebraic         := 7.4e+03
    }, $)::MATERIAL;

    description "Define a MATERIAL object with inputs: name of the material "
      "<name>, elastic modulus <elastic_modulus> (default = 210.0E9 (Pa)), "
      "Poisson's ratio <poisson_ratio> (default = 0.3 (-)), shear modulus "
      "<shear_modulus> (default = E/(2*(1+nu))), density <density> (default = "
      "7.4E3 (kg/m^3)).";

    return table({
      "type"            = MATERIAL,
      "name"            = name,
      "elastic_modulus" = elastic_modulus,
      "poisson_ratio"   = poisson_ratio,
      "shear_modulus"   = shear_modulus,
      "density"         = density
    });
  end proc: # MakeMaterial

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsMaterial := proc(
    obj::anything,
    $)::boolean;

    description "Check if the variable <var> is of MATERIAL type.";

    return type(obj, table) and evalb(var["type"] = MATERIAL);
  end proc: # IsMaterial

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeNode::static := proc({
    name::string,
    constrained_dof::list,
    coords::VECTOR,
    RF::FRAME := ground,
    }, $)::NODE;

    description "Define a MATERIAL object with inputs: name of the material "
      "<name>, elastic modulus <elastic_modulus> (default = 210.0E9 (Pa)), "
      "Poisson's ratio <poisson_ratio> (default = 0.3 (-)), shear modulus "
      "<shear_modulus> (default = E/(2*(1+nu))), density <density> (default = "
      "7.4E3 (kg/m^3)).";

    return table({
      "type"  = NODE,
      "name"  = name,
      "frame" = RF,
    });
  end proc: # MakeNode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsNode::static := proc(
    var::anything,
    $)::boolean;

    description "Check if the variable <var> is of NODE type.";

    return type(var, table) and evalb(var["type"] = NODE);
  end proc: # IsNode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/TrussMe_FEM/Affine.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # TrussMe_FEM

# That's all folks!
