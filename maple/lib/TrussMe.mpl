# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#           _____                   __  __        _____ _____ __  __          #
#          |_   _| __ _   _ ___ ___|  \/  | ___  |  ___| ____|  \/  |         #
#            | || '__| | | / __/ __| |\/| |/ _ \ | |_  |  _| | |\/| |         #
#            | || |  | |_| \__ \__ \ |  | |  __/ |  _| | |___| |  | |         #
#            |_||_|   \__,_|___/___/_|  |_|\___| |_|   |_____|_|  |_|         #
#        A package for mixed symbolic-numerical FEM structural analysis       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

TrussMe_FEM := module()

  description "A package for mixed symbolic-numerical FEM structural analysis.";

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  local m_IdLength     := 5;
  local m_VerboseMode  := false;
  local m_WarningMode  := true;
  local m_TimeLimit    := 2.0;
  local m_NodeColor    := "MediumSeaGreen";
  local m_SupportColor := "DarkOrange";
  local m_ElementColor := "SteelBlue";
  local m_ShellColor   := "MediumSlateBlue";
  local m_ForceColor   := "Firebrick";
  local m_MomentColor  := "Indigo";
  local m_NodeToken    := solidsphere;
  local m_SupportToken := solidbox;

  local m_CodegenOptions;
  local m_WorkingDirectory;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info := proc()

    description "Print module information.";

    printf(
      "+-------------------------------------------------------------------------------+\n"
      "| 'TrussMe[FEM]' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023 |\n"
      "| Current version authors:                                                      |\n"
      "|   Matteo Larcher and Davide Stocco.                                           |\n"
      "+-------------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad := proc()

    description "Module load procedure.";

    local lib_base_path, i, types;

    lib_base_path := NULL;
    for i in [libname] do
      if evalb(StringTools:-Search("TrussMe_FEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if evalb(lib_base_path = NULL) then
      error("cannot find 'TrussMe[FEM]' library in the 'TrussMe_FEM' toolbox "
        "folder.");
    end if;

    # Define module types
    TypeTools:-AddType('FRAME',      TrussMe_FEM:-IsFRAME);      protect('FRAME');
    TypeTools:-AddType('VECTOR',     TrussMe_FEM:-IsVECTOR);     protect('VECTOR');
    TypeTools:-AddType('POINT',      TrussMe_FEM:-IsPOINT);      protect('POINT');
    TypeTools:-AddType('MATERIAL',   TrussMe_FEM:-IsMATERIAL);   protect('MATERIAL');
    TypeTools:-AddType('DOFS',       TrussMe_FEM:-IsDOFS);       protect('DOFS');
    TypeTools:-AddType('NODE',       TrussMe_FEM:-IsNODE);       protect('NODE');
    TypeTools:-AddType('NODES',      TrussMe_FEM:-IsNODES);      protect('NODES');
    TypeTools:-AddType('SUPPORT',    TrussMe_FEM:-IsSUPPORT);    protect('SUPPORT');
    TypeTools:-AddType('SUPPORTS',   TrussMe_FEM:-IsSUPPORTS);   protect('SUPPORTS');
    TypeTools:-AddType('STIFFNESS',  TrussMe_FEM:-IsSTIFFNESS);  protect('STIFFNESS');
    TypeTools:-AddType('ELEMENT',    TrussMe_FEM:-IsELEMENT);    protect('ELEMENT');
    TypeTools:-AddType('ELEMENTS',   TrussMe_FEM:-IsELEMENTS);   protect('ELEMENTS');
    TypeTools:-AddType('STRUCTURE',  TrussMe_FEM:-IsSTRUCTURE);  protect('STRUCTURE');
    TypeTools:-AddType('COMPONENTS', TrussMe_FEM:-IsCOMPONENTS); protect('COMPONENTS');
    TypeTools:-AddType('LOAD',       TrussMe_FEM:-IsLOAD);       protect('LOAD');
    TypeTools:-AddType('LOADS',      TrussMe_FEM:-IsLOADS);      protect('LOADS');
    TypeTools:-AddType('FEM',        TrussMe_FEM:-IsFEM);        protect('FEM');

    # Codegen options
    TrussMe_FEM:-m_CodegenOptions := table([
      optimize          = true,
      digits            = 30,
      deducereturn      = false,
      coercetypes       = false,
      deducetypes       = false,
      reduceanalysis    = true,
      defaulttype       = numeric,
      functionprecision = double
    ]);
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload := proc()

    description "Module unload procedure.";

    unprotect('FRAME');      try TypeTools:-RemoveType('FRAME') catch: end try;
    unprotect('VECTOR');     try TypeTools:-RemoveType('VECTOR') catch: end try;
    unprotect('POINT');      try TypeTools:-RemoveType('POINT') catch: end try;
    unprotect('MATERIAL');   try TypeTools:-RemoveType('MATERIAL') catch: end try;
    unprotect('DOFS');       try TypeTools:-RemoveType('DOFS') catch: end try;
    unprotect('NODE');       try TypeTools:-RemoveType('NODE') catch: end try;
    unprotect('NODES');      try TypeTools:-RemoveType('NODES') catch: end try;
    unprotect('SUPPORT');    try TypeTools:-RemoveType('SUPPORT') catch: end try;
    unprotect('SUPPORTS');   try TypeTools:-RemoveType('SUPPORTS') catch: end try;
    unprotect('STIFFNESS');  try TypeTools:-RemoveType('STIFFNESS') catch: end try;
    unprotect('ELEMENT');    try TypeTools:-RemoveType('ELEMENT') catch: end try;
    unprotect('ELEMENTS');   try TypeTools:-RemoveType('ELEMENTS') catch: end try;
    unprotect('STRUCTURE');  try TypeTools:-RemoveType('STRUCTURE') catch: end try;
    unprotect('COMPONENTS'); try TypeTools:-RemoveType('COMPONENTS') catch: end try;
    unprotect('LOAD');       try TypeTools:-RemoveType('LOAD') catch: end try;
    unprotect('LOADS');      try TypeTools:-RemoveType('LOADS') catch: end try;
    unprotect('FEM');        try TypeTools:-RemoveType('FEM') catch: end try;

    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetModuleOptions := proc(
    {
    VerboseMode::{boolean, nothing}   := NULL,
    WarningMode::{boolean, nothing}   := NULL,
    TimeLimit::{nonnegative, nothing} := NULL,
    IdLength::{positive, nothing}     := NULL,
    NodeColor::{string, nothing}      := NULL,
    SupportColor::{string, nothing}   := NULL,
    ElementColor::{string, nothing}   := NULL,
    ShellColor::{string, nothing}     := NULL,
    ForceColor::{string, nothing}     := NULL,
    MomentColor::{string, nothing}    := NULL,
    NodeToken::{string, nothing}      := NULL,
    SupportToken::{string, nothing}   := NULL
    }, $)

    description "Set the module options: verbose mode <VerboseMode>, warning "
      "mode <WarningMode>, time limit <TimeLimit>, id code length <IdLength>, "
      "node color <NodeColor>, support color <SupportColor>, element color "
      "<ElementColor>, shell (2+ nodes) color <ShellColor>, force color "
      "<ForceColor>, moment color <MomentColor>, node token <NodeToken>, "
      "support token <SupportToken>.";

    if evalb(VerboseMode <> NULL) then
      TrussMe_FEM:-m_VerboseMode := VerboseMode;
    end if;

    if evalb(WarningMode <> NULL) then
      TrussMe_FEM:-m_WarningMode := WarningMode;
    end if;

    if evalb(TimeLimit <> NULL) then
      TrussMe_FEM:-m_TimeLimit := TimeLimit;
    end if;

    if evalb(IdLength <> NULL) then
      TrussMe_FEM:-m_IdLength := IdLength;
    end if;

    if evalb(NodeColor <> NULL) then
      TrussMe_FEM:-m_NodeColor := NodeColor;
    end if;

    if evalb(SupportColor <> NULL) then
      TrussMe_FEM:-m_SupportColor := SupportColor;
    end if;

    if evalb(ElementColor <> NULL) then
      TrussMe_FEM:-m_ElementColor := ElementColor;
    end if;

    if evalb(ShellColor <> NULL) then
      TrussMe_FEM:-m_ShellColor := ShellColor;
    end if;

    if evalb(ForceColor <> NULL) then
      TrussMe_FEM:-m_ForceColor := ForceColor;
    end if;

    if evalb(MomentColor <> NULL) then
      TrussMe_FEM:-m_MomentColor := MomentColor;
    end if;

    if evalb(NodeToken <> NULL) then
      TrussMe_FEM:-m_NodeToken := NodeToken;
    end if;

    if evalb(SupportToken <> NULL) then
      TrussMe_FEM:-m_SupportToken := SupportToken;
    end if;

    return NULL;
  end proc: # SetModuleOptions

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetObjByName := proc(
    objs::list(anything),
    name::string,
    $)::anything;

    description "Get object which field name is <name> from a list or set of "
      "objects <objs>.";

    local out, pos, i;

    out := NULL;
    for i from 1 to nops(objs) do
      if evalb(objs[i]["name"] = name) then
        out := eval(objs[i]); # Do not remove eval: eval(table)
        pos := i;
        break;
      end if;
    end do;

    if (_nresults = 2) then
      return eval(out), pos;
    else
      return eval(out);
    end if;
  end proc: # GetObjByName

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetObjById := proc(
    objs::list(anything),
    id_fld::string,
    {
    position::boolean := false
    }, $)::anything;

    description "Get object which field 'id' is equal to <id_fld>, between "
      "the objects <objs>.";

    local out, pos, i;

    out := NULL;
    for i from 1 to nops(objs) do
      if evalb(objs[i]["id"] = id_fld) then
        out := `if`(position, i, eval(objs[i])); # Do not remove eval: eval(table)
        break;
      end if;
    end do;
    return eval(out);
  end proc: # GetObjByName

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetObjsByType := proc(
    objs::list(anything),
    type_fld::{symbol, list(symbol), set(symbol)},
    {
    position::boolean := false
    }, $)::list(anything);

    description "Get objects which field 'type' is equal to <type_fld>, between "
      "the objects <objs>.";

    local type_set, out, i;

    if type(type_fld, symbol) then
      type_set := {type_fld};
    elif type(type_fld, list) then
      type_set := convert(type_fld, set);
    else
      type_set := type_fld;
    end if;

    out := [];
    for i from 1 to nops(objs) do
      if type(objs[i], type_set) then
        out := [op(out), `if`(position, i, eval(objs[i]))]; # Do not remove eval: eval(table)
      end if;
    end do;
    return eval(out);
  end proc: # GetObjsByType

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Simplify := proc(
    var::anything,
    opt::anything := NULL,
    $)::anything;

    description "Try to simplify an algebraic expression <var> with optional "
      "simplification options <opt>. The simplification is performed within "
      "the internal or indexed time limit.";

    try
      return timelimit(
        `if`(procname::indexed, op(procname), TrussMe_FEM:-m_TimeLimit),
        (simplify(var, opt) assuming real)
      );
    catch "time expired":
      if TrussMe_FEM:-m_WarningMode then
        WARNING("exceeded time limit, raw solutions is returned.");
      end if;
      return var;
    catch "division by zero":
      error("division by zero detected.");
      return var;
    catch:
      error("something went wrong, last exception '%1'.", lastexception);
      return var;
    end try:
  end proc: # Simplify

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Norm2 := proc(
    x::{list, Vector},
    $)::algebraic;

    description "Compute the Euclidean norm of list or vector <x>.";

    return sqrt(add(i, i in x^~2));
  end proc: # Norm2

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateId := proc({
      size::positive := TrussMe_FEM:-m_IdLength,
      opts::symbol   := 'alnum'
    }, $)::string;

    return StringTools:-Random(size, opts);
  end proc: # GenerateId

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Spy := proc(
    A::Matrix,
    $)::anything;

    description "Plot of non-zero values of matrix <A>.";

    return plots:-sparsematrixplot(A, 'matrixview');
  end proc: # Spy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/FEM/Affine.mpl"
$include "./lib/FEM/Codegen.mpl"
$include "./lib/FEM/Differentiation.mpl"
$include "./lib/FEM/Load.mpl"
$include "./lib/FEM/Material.mpl"
$include "./lib/FEM/Plot.mpl"
$include "./lib/FEM/Structure.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # TrussMe_FEM

# That's all folks!
