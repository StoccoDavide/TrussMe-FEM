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
    TypeTools:-AddType('FRAME',    TrussMe_FEM:-IsFRAME);    protect('FRAME');
    TypeTools:-AddType('VECTOR',   TrussMe_FEM:-IsVECTOR);   protect('VECTOR');
    TypeTools:-AddType('POINT',    TrussMe_FEM:-IsPOINT);    protect('POINT');
    TypeTools:-AddType('EARTH',    TrussMe_FEM:-IsEARTH);    protect('EARTH');
    TypeTools:-AddType('NODE',     TrussMe_FEM:-IsNODE);     protect('NODE');
    TypeTools:-AddType('SUPPORT',  TrussMe_FEM:-IsSUPPORT);  protect('SUPPORT');
    TypeTools:-AddType('DOFS',     TrussMe_FEM:-IsDOFS);     protect('DOFS');
    TypeTools:-AddType('ELEMENT',  TrussMe_FEM:-IsELEMENT);  protect('ELEMENT');
    TypeTools:-AddType('FORCE',    TrussMe_FEM:-IsFORCE);    protect('FORCE');
    TypeTools:-AddType('QFORCE',   TrussMe_FEM:-IsQFORCE);   protect('QFORCE');
    TypeTools:-AddType('MOMENT',   TrussMe_FEM:-IsMOMENT);   protect('MOMENT');
    TypeTools:-AddType('QMOMENT',  TrussMe_FEM:-IsQMOMENT);  protect('QMOMENT');
    TypeTools:-AddType('MATERIAL', TrussMe_FEM:-IsMATERIAL); protect('MATERIAL');
    #TypeTools:-AddType('STRUCTURE', TrussMe_FEM:-IsStructure); protect('STRUCTURE')
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()

    description "TrussMe_FEM' module unload procedure.";

    unprotect(
      'FRAME',
      'VECTOR',
      'POINT',
      'EARTH',
      'NODE',
      'SUPPORT',
      'DOFS',
      'ELEMENT',
      'FORCE',
      'QFORCE',
      'MOMENT',
      'QMOMENT',
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
    if _self:-WarningMode then
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

  export IsEARTH::static := proc(
    var::anything,
    $)::boolean;

    description "Check if the variable <var> is of EARTH type.";

    return type(var, table) and evalb(var["type"] = EARTH);
  end proc: # IsEARTH

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Simplify::static := proc(
    _self::TrussMe_FEM,
    var::anything,
    opt::anything := NULL,
    $)::anything;

    description "Try to simplify an algebraic expression <var> with optional "
      "simplification options <opt>. The simplification is performed within "
      "the internal or indexed time limit of <TimeLimit> seconds.";

    try
      return timelimit(
        `if`(procname::indexed, op(procname), TrussMe_FEM:-m_TimeLimit),
        simplify(var, opt)
      );
    catch "time expired":
      if _self:-WarningMode then
        WARNING("time limit of exceeded, raw solutions is returned.");
      end if;
      return var;
    catch "division by zero":
      error("division by zero detected.");
      return var;
    catch:
      error("something went wrong.");
    end try:
  end proc: # Simplify

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateId::static := proc({
      size::positive := 5,
      opts::symbol   := 'alnum'
    }, $)::string;

    return StringTools:-Random(size, opts);
  end proc: # GenerateId

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/TrussMe_FEM/Affine.mpl"
$include "./lib/TrussMe_FEM/Element.mpl"
$include "./lib/TrussMe_FEM/Loads.mpl"
$include "./lib/TrussMe_FEM/Material.mpl"
$include "./lib/TrussMe_FEM/Node.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # TrussMe_FEM

# That's all folks!
