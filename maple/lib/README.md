# MAPLE API

```
# A package for mixed symbolic-numerical FEM structural analysis.
module :

    # Print module information.
    Info( )

    # Module load procedure.
    ModuleLoad( )

    # Module unload procedure.
    ModuleUnload( )

    # Set the module options: verbose mode <VerboseMode>, warning mode
    # <WarningMode>, time limit <TimeLimit>, id code length <IdLength>, node
    # color <NodeColor>, support color <SupportColor>, element color
    # <ElementColor>, shell (2+ nodes) color <ShellColor>, force color
    # <ForceColor>, moment color <MomentColor>, node token <NodeToken>, support
    # token <SupportToken>.
    SetModuleOptions( { ElementColor::{nothing, string} := NULL,
                        ForceColor::{nothing, string} := NULL,
                        IdLength::{nothing, positive} := NULL,
                        MomentColor::{nothing, string} := NULL,
                        NodeColor::{nothing, string} := NULL,
                        NodeToken::{nothing, string} := NULL,
                        ShellColor::{nothing, string} := NULL,
                        SupportColor::{nothing, string} := NULL,
                        SupportToken::{nothing, string} := NULL,
                        TimeLimit::{nothing, nonnegative} := NULL,
                        VerboseMode::{boolean, nothing} := NULL,
                        WarningMode::{boolean, nothing} := NULL }, $ )

    # Set gravity vector with [x, y, z]^T components of <vec>.
    SetGravity( g::VECTOR, $ )

    # Get gravity vector.
    GetGravity( $ )

    # Get object which field name is <name> from a list or set of objects
    # <objs>.
    GetObjByName( objs::list(anything), name::string, $ ) :: anything

    # Get object which field 'id' is equal to <id_fld>, between the objects
    # <objs>.
    GetObjById( objs::list(anything), id_fld::string,
                { position::boolean := false }, $ ) :: anything

    # Get objects which field 'type' is equal to <type_fld>, between the
    # objects <objs>.
    GetObjsByType( objs::list(anything),
                   type_fld::{symbol, list(symbol), set(symbol)},
                   { position::boolean := false }, $ ) :: list(anything)

    # Try to simplify an algebraic expression <var> with optional
    # simplification options <opt>. The simplification is performed within the
    # internal or indexed time limit.
    Simplify( var::anything, opt::anything := NULL, $ ) :: anything

    # Compute the Euclidean norm of list or vector <x>.
    Norm2( x::{Vector, list}, $ ) :: algebraic

    GenerateId( { opts::symbol := 'alnum',
                  size::positive := TrussMe_FEM:-m_IdLength }, $ ) :: string

    # Plot of non-zero values of matrix <A>.
    Spy( A::Matrix, $ ) :: anything

    # Check if the variable <var> is of FRAME type.
    IsFRAME( var::anything, $ ) :: boolean

    # Generate a reference frame matrix from three points or vectors <p_1>,
    # <p_2> and vector <vec> ortogonal to XY-plane <vec>
    GenerateFrame( p_1::{POINT, VECTOR, Vector, list},
                   p_2::{POINT, VECTOR, Vector, list},
                   vec::{VECTOR, Vector, list}, $ ) :: FRAME

    # Generate a generic reference frame matrix from a string label <label>.
    # Optional arguments <e>, <p>, <x>, <y>, <z> and <s> are used to customize
    # the output, e.g., 'exx := e||s||x||x'.
    GenerateGenericFrame( label::string := "",
                          { e := "e", p := "p", s := "__", x := "x", y := "y",
                            z := "z" }, $ ) :: FRAME

    # Inverse affine transformation matrix <RF>.
    InverseFrame( RF::FRAME, $ ) :: FRAME

    # Affine transformation matrix corresponding to a translation <x, y, z>.
    Translate( x::algebraic, y::algebraic, z::algebraic, $ ) :: FRAME

    # Extract the translation vector of the reference frame <RF>.
    Translation( RF::FRAME, $ ) :: Vector

    # Affine transformation matrix corresponding to the rotation <angle> around
    # the given <axis>.
    Rotate( axis::{string, symbol}, angle::algebraic, $ ) :: FRAME

    # Extract the rotation matrix of the reference frame <RF>.
    Rotation( RF::FRAME, $ ) :: Matrix

    # Check if the variable <var> is of VECTOR type.
    IsVECTOR( var::anything, $ ) :: boolean

    # Check if the variable <var> is of POINT type.
    IsPOINT( var::anything, $ ) :: boolean

    # Extract the origin point of the reference frame <RF>.
    Origin( RF::FRAME, $ ) :: POINT

    # Extract the x-axis component of the vector or point <x>.
    CompX( x::{POINT, VECTOR}, $ ) :: algebraic

    # Extract the y-axis component of the vector or point <x>.
    CompY( x::{POINT, VECTOR}, $ ) :: algebraic

    # Extract the z-axis component of the vector or point <x>.
    CompZ( x::{POINT, VECTOR}, $ ) :: algebraic

    # Extract the x, y and z-axis components of the vector or point <x>.
    CompXYZ( x::{POINT, VECTOR}, $ ) :: algebraic

    # Extract the x-axis unit vector of the reference frame <RF>.
    UvecX( RF::FRAME := Matrix(4,shape = identity), $ ) :: VECTOR

    # Extract the y-axis unit vector of the reference frame <RF>.
    UvecY( RF::FRAME, $ ) :: VECTOR

    # Extract the z-axis unit vector of the reference frame <RF>.
    UvecZ( RF::FRAME := Matrix(4,shape = identity), $ ) :: VECTOR

    # Extract the x, y and z-axis unit vectors of the reference frame <RF>.
    UvecXYZ( RF::FRAME := Matrix(4,shape = identity), $ ) :: VECTOR

    # Project the vector or point <x>, from reference frame <RF_ini> to
    # reference frame <RF_end>.
    Project( x::{POINT, VECTOR}, RF_ini::FRAME, RF_end::FRAME, $ )
           :: {POINT, VECTOR}

    # Check if the variable <var> is of MATERIAL type.
    IsMATERIAL( var::anything, $ ) :: boolean

    # Define material with inputs: name <name>, elastic modulus (E)
    # <elastic_modulus>, Poisson's ratio (nu) <poisson_ratio>, shear modulus
    # <shear_modulus> (default = E/(2*(1+nu))), and density <density>.
    MakeMaterial( { density::algebraic := 0, elastic_modulus::algebraic := 0,
                    name::string := "Undefined",
                    poisson_ratio::algebraic := 0,
                    shear_modulus::algebraic :=
                      elastic_modulus/(2+2*poisson_ratio) },
                  $ ) :: MATERIAL

    # Get default steel material with 'CarbonSteel' name, elastic modulus E =
    # 210.0e+09 (Pa), Poisson's ratio nu = 0.3 (-), shear modulusE/(2*(1+nu),
    # density rho = 7850.0 (kg/m^3).
    MakeCarbonSteel( $ )

    # Get default steel material with 'InoxSteel' name, elastic modulus E =
    # 200.0e+09 (Pa), Poisson's ratio nu = 0.3 (-), shear modulusE/(2*(1+nu),
    # density rho = 8000.0 (kg/m^3).
    MakeInoxSteel( $ )

    # Get default titanium material with 'Titanium' name, elastic modulus E =
    # 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 4500.0 (kg/m^3).
    MakeTitanium( $ )

    # Get default copper material with 'Copper' name, elastic modulus E =
    # 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 4500.0 (kg/m^3).
    MakeCopper( $ )

    # Get default brass material with 'Brass' name, elastic modulus E =
    # 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 4500.0 (kg/m^3).
    MakeBrass( $ )

    # Get default bronze material with 'Bronze' name, elastic modulus E =
    # 110.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 4500.0 (kg/m^3).
    MakeBronze( $ )

    # Get default lead material with 'Lead' name, elastic modulus E = 110.0e+09
    # (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu), density
    # rho = 4500.0 (kg/m^3).
    MakeLead( $ )

    # Get default zinc material with 'Zinc' name, elastic modulus E = 110.0e+09
    # (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu), density
    # rho = 4500.0 (kg/m^3).
    MakeZinc( $ )

    # Get default magnesium material with 'Magnesium' name, elastic modulus E =
    # 45.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 1800.0 (kg/m^3).
    MakeMagnesium( $ )

    # Get default alluminium material with 'Alluminium' name, elastic modulus E
    # = 69.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 8000.0 (kg/m^3).
    MakeAlluminium( $ )

    # Get default avional material with 'Avional' name, elastic modulus E =
    # 70.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 2690.0 (kg/m^3).
    MakeAvional( $ )

    # Get default peraluman material with 'Peraluman' name, elastic modulus E =
    # 70.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 2700.0 (kg/m^3).
    MakePeraluman( $ )

    # Get default anticorodal material with 'Anticorodal' name, elastic modulus
    # E = 69.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear
    # modulusE/(2*(1+nu), density rho = 2700.0 (kg/m^3).
    MakeAnticorodal( $ )

    # Get default carpental material with 'Carpental' name, elastic modulus E =
    # 72.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 2780.0 (kg/m^3).
    MakeCarpental( $ )

    # Get default ergal material with 'Ergal' name, elastic modulus E =
    # 72.0e+09 (Pa), Poisson's ratio nu = 0.33 (-), shear modulusE/(2*(1+nu),
    # density rho = 2780.0 (kg/m^3).
    MakeErgal( $ )

    # Check if the variable <var> is of NODE type.
    IsNODE( var::anything, $ ) :: boolean

    # Check if the variable <var> is of NODES type.
    IsNODES( var::anything, $ ) :: boolean

    # Check if the variable <var> is of SUPPORT type.
    IsSUPPORT( var::anything, $ ) :: boolean

    # Check if the variable <var> is of SUPPORTS type.
    IsSUPPORTS( var::anything, $ ) :: boolean

    # Check if the variable <var> is of DOFS type.
    IsDOFS( var::anything, $ ) :: boolean

    # Create a node with name <name> at coordinates <coordinates> in the
    # reference frame <frame>. The constraints on dofs are specified by <dofs>,
    # where 1 means free and 0 means that the dof is constrained in the
    # direction given from (global to local) <frame>.
    MakeNode( name::string, coordinates::{POINT, VECTOR, Vector},
              { displacements::Vector(algebraic) := <0, 0, 0, 0, 0, 0>,
                dofs::DOFS := <1, 1, 1, 1, 1, 1>,
                frame::FRAME := Matrix(4,shape = identity) }, $ ) :: NODE

    # Create a node with name <name> at coordinates <coordinates> in the
    # reference frame <frame>. The constraints on dofs are specified by <dofs>,
    # where 1 means free and 0 means that the dof is constrained in the
    # direction given from (global to local) <frame>. The node is also
    # connected to a compliant spring element with traslational stiffness <K>
    # and torsional stiffness <T>.
    MakeCompliantNode( name::string, coordinates::{POINT, VECTOR, Vector},
                       { K::{algebraic, list(algebraic)} := 0,
                         T::{algebraic, list(algebraic)} := 0,
                         displacements::Vector(algebraic) := <0, 0, 0, 0, 0,
                           0>,
                         dofs::DOFS := <1, 1, 1, 1, 1, 1>,
                         frame::FRAME := Matrix(4,shape = identity) }, $ )
                     :: NODE

    # Check if the variable <var> is of STIFFNESS type.
    IsSTIFFNESS( var::anything, $ ) :: boolean

    # Get the stiffness matrix of a spring given the translational stiffnesses
    # <K> or <K_x, K_y, K_z>, and the torsional stiffnesses <T> or <T_x, T_y,
    # T_z>.
    GetSpringStiffness( K::{algebraic, list(algebraic)},
                        T::{algebraic, list(algebraic)}, $ ) :: STIFFNESS

    # Get the stiffness matrix of a rod (only axial-stiffness) given the
    # cross-section area <A>, the Young's modulus <E>, the shear modulus <G>,
    # the length <L>, and the cross-section inertia <I_x>, <I_y> and <I_z>.
    GetRodStiffness( A::algebraic, E::algebraic, G::algebraic, L::algebraic,
                     I_x::algebraic, I_y::algebraic, I_z::algebraic, $ )
                   :: STIFFNESS

    # Get the stiffness matrix of a lean beam given the cross-section area <A>,
    # the Young's modulus <E>, the shear modulus <G>, the length <L>, and the
    # inertia of the cross-section <I_x>, <I_y> and <I_z>.
    GetBeamStiffness( A::algebraic, E::algebraic, G::algebraic, L::algebraic,
                      I_x::algebraic, I_y::algebraic, I_z::algebraic, $ )
                    :: STIFFNESS

    # Get the stiffness matrix of a Timoshenko's (thick) beam given the
    # cross-section area <A>, the Young's modulus <E>, the shear modulus <G>,
    # the length <L>, and the inertia of the cross-section <I_x>, <I_y> and
    # <I_z>.
    GetTimoshenkoBeamStiffness( A::algebraic, E::algebraic, G::algebraic,
                                L::algebraic, I_x::algebraic, I_y::algebraic,
                                I_z::algebraic, $ ) :: STIFFNESS

    # Get the output displacements of a solved and stored FEMstructure <fem>.
    GetOutputDisplacements( fem::_FEM, $ ) :: Matrix

    # Get the output reactions of a solved and stored FEMstructure <fem>.
    GetOutputReactions( fem::_FEM, $ ) :: Matrix

    # Check if the variable <var> is of ELEMENT type.
    IsELEMENT( var::anything, $ ) :: boolean

    # Check if the variable <var> is of ELEMENTS type.
    IsELEMENTS( var::anything, $ ) :: boolean

    # Make an element with name <name> on reference frame <frame>, connecting
    # the dofs <Ni_dofs> on i-th node i <[Ni, Ni_dofs]>, connecting with
    # stiffness <K>.
    MakeElement( name::string,
                 { frame::FRAME := TrussMe_FEM:-GenerateGenericFrame(name) } )
               :: ELEMENT

    # Make a spring element with name <name> on reference frame <frame>,
    # connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting the
    # dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with translational stiffnesses
    # <K> or <K_x, K_y, K_z>, and torsional stiffnesses <T> or <T_x, T_y, T_z>.
    # Optional nodes distance <distance> can be specified.
    MakeSpring( name::string, N1::{NODE, list({DOFS, NODE})},
                N2::{NODE, list({DOFS, NODE})},
                { K::{algebraic, list(algebraic)} := 0,
                  T::{algebraic, list(algebraic)} := 0,
                  distance::algebraic := 0,
                  frame::FRAME := TrussMe_FEM:-GenerateGenericFrame(name) },
                $ ) :: ELEMENT

    # Make a rod element with name <name> on reference frame <frame>,
    # connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting the
    # dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with material <material>,
    # cross-section area <area>and inertia <inertia>. Optional nodes distance
    # <distance> can be specified.
    MakeRod( name::string,
             N1::{NODE, Vector({DOFS, NODE}), list({DOFS, NODE})},
             N2::{NODE, Vector({DOFS, NODE}), list({DOFS, NODE})},
             { area::algebraic := 0, distance::algebraic := -1,
               frame::FRAME := TrussMe_FEM:-GenerateGenericFrame(name),
               inertia::list(algebraic) := [0, 0, 0],
               material::MATERIAL := TrussMe_FEM:-MakeCarbonSteel() }, $ )
           :: ELEMENT

    # Make a beam element with name <name> on reference frame <frame>,
    # connecting the dofs <N1_dofs> on node 1 <[N1, N1_dofs]>, connecting the
    # dofs <N2_dofs> on node 2 <[N2, N2_dofs]> with material <material> and
    # cross-section area <area> and inertia <inertia>. Optional nodes distance
    # <distance> and Timoshenko's beam boolean <timoshenko> can be specified.
    MakeBeam( name::string,
              N1::{NODE, Vector({DOFS, NODE}), list({DOFS, NODE})},
              N2::{NODE, Vector({DOFS, NODE}), list({DOFS, NODE})},
              { area::algebraic := 0, distance::algebraic := -1,
                frame::FRAME := TrussMe_FEM:-GenerateGenericFrame(name),
                inertia::list(algebraic) := [0, 0, 0],
                material::MATERIAL := TrussMe_FEM:-MakeCarbonSteel(),
                timoshenko::boolean := false }, $ ) :: ELEMENT

    # Check if the variable <var> is of STRUCTURE type.
    IsSTRUCTURE( var::anything, $ ) :: boolean

    # Get nodal dofs of nodes <nodes>.
    GetNodalDofs( nodes::NODES, $ ) :: Vector

    # Get the vector of nodal displacements of nodes <nodes>.
    GetNodalDisplacements( nodes::NODES, $ ) :: Vector

    # Compute the transformation matrix of the nodes <nodes>.
    StiffnessTransformation( nodes::NODES, $ ) :: Matrix

    # Compute the global stiffness matrix of the nodes <nodes> and elements
    # <elements>.
    GlobalStiffness( nodes::NODES, elements::ELEMENTS, $ ) :: Matrix

    # Compute the global stiffness matrix of the nodes <nodes> and elements
    # <elements>.
    GlobalStiffnessPrime( nodes::NODES, elements::ELEMENTS, $ ) :: Matrix

    # Check if the variable <var> is of FEM type.
    IsFEM( var::anything, $ ) :: boolean

    # Split the FEM structure <fem> in free and constrained dofs.
    SplitFEM( fem::_FEM, $ )

    # Recast the system <fem> to avoid singularities.
    RecastFEM( fem::_FEM, $ )

    # Generate a FEM structure from the nodes <nodes>, elements <elements>, and
    # loads <loads>. If <tryhard> is enabled, the stiffness matrix will be
    # recasted to avoid singularities.
    GenerateFEM( nodes::NODES, elements::ELEMENTS, loads::LOADS,
                 { tryhard::boolean := false }, $ ) :: _FEM

    # Solve the FEM structure <fem> and optionally use LAST LU decompostion
    # <use_LAST> and veil the expressions <use_LEM> with signature mode
    # <use_SIG> and label <label>. Factorization method <factorization> can be
    # choosen between 'LU', fraction-free 'FFLU', 'QR', and Gauss-Jordan 'GJ'.
    SolveFEM( fem::_FEM,
              { factorization::string := "LU", label::string := "V",
                use_LAST::boolean := false, use_LEM::boolean := true,
                use_SIG::boolean := true }, $ )

    # Store the FEM structure <fem> in the nodes <nodes>.
    StoreFEM( fem::_FEM, $ )

    # Check if the variable <var> is of LOAD type.
    IsLOAD( var::anything, $ ) :: boolean

    # Check if the variable <var> is of COMPONENTS type.
    IsCOMPONENTS( var::anything, $ ) :: boolean

    # Create a load with name <name> acting on the node (or its id) <node> on
    # the frame <frame> with components <components>.
    MakeLoad( name::string, node::NODE, components::COMPONENTS,
              { frame::{FRAME, string} := node["id"] }, $ ) :: LOAD

    # Check if the variable <var> is of LOADS type.
    IsLOADS( var::anything, $ ) :: boolean

    # Get the vector of nodal loads <loads> of nodes <nodes>.
    GetNodalLoads( nodes::NODES, loads::LOADS, $ ) :: Vector

    # Generate a file named <fname> with the content <str>.
    GenerateFile( fname::string, str::string, $ )

    # Procedure that creates an empty file <fname>.
    ClearFile( fname::string, $ )

    # Set options for code generation optimization.
    SetCodegenOptions( { coercetypes::boolean :=
                           TrussMe_FEM:-m_CodegenOptions[parse("coercetypes")],
                         deducereturn::boolean :=
                           TrussMe_FEM:-m_CodegenOptions[parse("deducereturn")],
                         deducetypes::boolean :=
                           TrussMe_FEM:-m_CodegenOptions[parse("deducetypes")],
                         defaulttype::symbol :=
                           TrussMe_FEM:-m_CodegenOptions[parse("defaulttype")],
                         digits::posint :=
                           TrussMe_FEM:-m_CodegenOptions[parse("digits")],
                         functionprecision::symbol :=
                           TrussMe_FEM:-m_CodegenOptions[parse("functionprecision")],
                         optimize::boolean :=
                           TrussMe_FEM:-m_CodegenOptions[parse("optimize")],
                         reduceanalysis::boolean :=
                           TrussMe_FEM:-m_CodegenOptions[parse("reduceanalysis")]
                           },
                       $ )

    # Get options for code generation.
    GetCodegenOptions( field::string := "all", $ ) :: anything

    # Convert a list of expressions <expr_list> into Matlab code.
    TranslateToMatlab( expr_list::list, $ )

    # Apply indentation <ind> to string <str>.
    ApplyIndent( ind::string, str::string, $ ) :: string

    # Generate properties code from a list of data <data> and optional
    # indentation <indent>.
    GenerateProperties( data::list(symbol),
                        { indent::string := "  " }, $ ) :: string

    # Generate inputs code from a list of variables <vars> and optional
    # indentation <indent> and skip null inputs <skipnull>.
    GenerateInputs( vars::list(list({string, symbol})),
                    { indent::string := "  ", skipnull::boolean := true }, $ )
                  :: string

    # Extract elements for a n-dimensional function <func> with name <name>,
    # dimensions <dims>, and optional veiling label <label>, indentation
    # <indent> and skip null inputs <skipnull>.
    ExtractElements( name::string, func::{Array, Matrix, Vector, list},
                     dims::list(nonnegint),
                     { indent::string := "  ", label::string := "out",
                       skipnull::boolean := true }, $ ) :: list

    # Generate a function header for a function <name> with variables <vars>,
    # optional description <info> and indentation <indent> and skip class
    # object input <skipthis>.
    GenerateHeader( name::string, vars::list(list(symbol)),
                    { indent::string := "  ",
                      info::string := "No description provided.",
                      skipthis::boolean := false }, $ ) :: string

    # Generate code for elements <func> with optional indentation <indent>.
    GenerateElements( func::{Array, Matrix, Vector, list},
                      { indent::string := "  " }, $ ) :: string

    # Generate code for function body for a function <name> with dimensions
    # <dims>, optional header <header>, properties <properties>, inputs
    # <inputs>, veils <veils>, elements <elements>, indentation <indent>,
    # outputs <outputs> and vector type <typestr>.
    GenerateBody( name::string, dims::list(nonnegint),
                  { elements::string := "No elements",
                    header::string := "No header", indent::string := "  ",
                    inputs::string := "No inputs",
                    outputs::string := "No outputs",
                    properties::string := "No properties",
                    typestr::string := "zeros", veils::string := "No veils" },
                  $ ) :: string

    # Translate the vector <vec> with variables <vars> into a Matlab function
    # named <name> and return it as a string. The optional arguments and class
    # properties <data>, function description <info>, veiling label <label>,
    # and indentation string <indent>.
    VectorToMatlab( name::string, vars::list(list(symbol)), vec::Vector,
                    { data::list(symbol) := [], indent::string := "  ",
                      info::string := "No info", label::string := "out",
                      skipnull::boolean := true }, $ ) :: string

    # Translate the vector <vec> with variables <vars> into a Matlab function
    # named <name> and return it as a string. The optional arguments and class
    # properties <data>, function description <info>, veiling label <label>,
    # and indentation string <indent>.
    MatrixToMatlab( name::string, vars::list(list(symbol)), mat::Matrix,
                    { data::list(symbol) := [], indent::string := "  ",
                      info::string := "No info", label::string := "out",
                      skipnull::boolean := true }, $ ) :: string

    # Get the substitutions to transform veils <vars> from (v[j])(f) to v_j.
    GetVeilSubs( veil::{Vector(algebraic), list(algebraic)}, $ )

    # Generate a constructor for a system named <name> system data <data>,
    # description <info> and indentation <indent>.
    GenerateConstructor( name::string,
                         { data::list(symbol = algebraic) := [],
                           indent::string := "  ",
                           info::string := "Class constructor." }, $ )
                       :: string

    # Generate a FEM system <fem> with name <name>, with optional data <data>,
    # description <info>, output label <label> and indentation <indent>.
    SystemToMatlab( name::string, fem::_FEM,
                    { data::list(symbol = algebraic) := [],
                      indent::string := "  ",
                      info::string := "No class description provided.",
                      label::string := "out", vars::list(symbol) := [] }, $ )
                  :: string

    # Generate Matlab code for the FEM system <fem> with name <name>, with
    # optional data <data>, description <info>, output label <label> and
    # indentation <indent>, variables <vars> and output path <path>.
    GenerateMatlabCode( name::string, fem::_FEM,
                        { data::list(symbol = algebraic) := [],
                          indent::string := "  ",
                          info::string := "No class description provided.",
                          label::string := "out", path::string := "./",
                          vars::list(symbol) := [] }, $ )

    # Return the color of the object <obj>.
    ObjectColor( obj::{ELEMENT, LOAD, NODE, SUPPORT}, $ ) :: string

    # Return the graph of the structure given a list of nodes <nodes> and
    # elements <elements>.
    StructureGraph( fem::_FEM,
                    { disp::boolean := true, id::boolean := false }, $ )
                  :: function

    # Draw a reference frame <frame> given a list of substitution data <data>,
    # an axes scaling factor <scaling>, and axes colors <colors>.
    DrawFrame( frame::FRAME,
               { colors::list(string) := ["Red", "Green", "Blue"],
                 data::{list(`=`), set(`=`)} := [], scaling::numeric := 1.0 },
               $ ) :: function

    # Plot the node (or support) at point <p> given a list or set of data for
    # substitution <data>, a display token <token> and a display color <color>.
    PlotNode( p::{Vector(algebraic), list(algebraic)},
              { color::string := TrussMe_FEM:-m_NodeColor,
                data::{list(`=`), set(`=`)} := [],
                token::symbol := TrussMe_FEM:-m_NodeToken }, $ ) :: function

    # Plot the element from point <p_1> and <p_2> given a list or set of data
    # for substitution <data> and a display color <color>.
    PlotElement( p_1::{Vector(algebraic), list(algebraic)},
                 p_2::{Vector(algebraic), list(algebraic)},
                 { color::string := TrussMe_FEM:-m_ElementColor,
                   data::{list(`=`), set(`=`)} := [] }, $ ) :: function

    # Plot the element from diplacements <d_1> and <d_2> given a list or set of
    # data for substitution <data> and a display color <color>.
    PlotDeformedElement( p_1::{Vector(algebraic), list(algebraic)},
                         p_2::{Vector(algebraic), list(algebraic)},
                         d_1::{Vector(algebraic), list(algebraic)},
                         d_2::{Vector(algebraic), list(algebraic)},
                         { color::string := TrussMe_FEM:-m_ElementColor,
                           data::{list(`=`), set(`=`)} := [],
                           scaling::nonnegative := 1.0 }, $ ) :: function

    # Plot the load arrow from point <p_1> and <p_2> given a list or set of
    # data for substitution <data> and a display color <color>.
    PlotLoad( p_1::{Vector(algebraic), list(algebraic)},
              p_2::{Vector(algebraic), list(algebraic)},
              { color::string := TrussMe_FEM:-m_LoadColor,
                data::{list(`=`), set(`=`)} := [],
                scaling::nonnegative := 1.0 }, $ ) :: function

    # Plot the undeformed <fem> structure given a list or set of substitution
    # data <data>, a frame scaling factor <frame_scaling> or <nodes_scaling,
    # elements_scaling>, and a loads scaling factor <scaling>.
    PlotStructure( fem::_FEM,
                   { data::{list(`=`), set(`=`)} := [],
                     frame_scaling::{numeric, list(numeric)} := 0.,
                     load_scaling::numeric := 1.0 }, $ )
                 :: {function, list(function)}

    # Plot the deformed <fem> structure given a list or set of substitution
    # data <data>, a frame scaling factor <frame_scaling> or <nodes_scaling,
    # elements_scaling>, a loads scaling factor <load_scaling>, and a
    # deformation magnification factor <deformation_scaling>.
    PlotDeformedStructure( fem::_FEM,
                           { data::{list(`=`), set(`=`)} := [],
                             deformation_scaling::numeric := 1.0,
                             frame_scaling::{numeric, list(numeric)} := 0.,
                             interpolate::boolean := true,
                             load_scaling::numeric := 1.0 }, $ )
                         :: {function, list(function)}
```
