# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                    ____          _                                          #
#                   / ___|___   __| | ___  __ _  ___ _ __                     #
#                  | |   / _ \ / _` |/ _ \/ _` |/ _ \ '_ \                    #
#                  | |__| (_) | (_| |  __/ (_| |  __/ | | |                   #
#                   \____\___/ \__,_|\___|\__, |\___|_| |_|                   #
#                                         |___/                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateFile := proc(
  fname::string,
  str::string,
  $)

  description "Generate a file named <fname> with the content <str>.";

  local file;

  file := FileTools:-Text:-Open(fname, create = true, overwrite = true);
  FileTools:-Text:-WriteString(file, str);
  FileTools:-Text:-Close(file);
  return NULL;
end proc: # GenerateFile

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ClearFile := proc(
  fname::string,
  $)

  description "Procedure that creates an empty file <fname>.";

  local file;

  file := FileTools:-Text:-Open(fname, create = true, overwrite = true);
  FileTools:-Text:-WriteString(file, "\n");
  FileTools:-Text:-Close(file);
  return NULL;
end proc: # ClearFile

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SetCodegenOptions := proc({
    optimize::boolean         := TrussMe:-FEM:-m_CodegenOptions[parse("optimize")],
    digits::posint            := TrussMe:-FEM:-m_CodegenOptions[parse("digits")],
    deducereturn::boolean     := TrussMe:-FEM:-m_CodegenOptions[parse("deducereturn")],
    coercetypes::boolean      := TrussMe:-FEM:-m_CodegenOptions[parse("coercetypes")],
    deducetypes::boolean      := TrussMe:-FEM:-m_CodegenOptions[parse("deducetypes")],
    reduceanalysis::boolean   := TrussMe:-FEM:-m_CodegenOptions[parse("reduceanalysis")],
    defaulttype::symbol       := TrussMe:-FEM:-m_CodegenOptions[parse("defaulttype")],
    functionprecision::symbol := TrussMe:-FEM:-m_CodegenOptions[parse("functionprecision")]
  }, $)

  description "Set options for code generation optimization.";

  TrussMe:-FEM:-m_CodegenOptions[parse("optimize")]          := optimize;
  TrussMe:-FEM:-m_CodegenOptions[parse("digits")]            := digits;
  TrussMe:-FEM:-m_CodegenOptions[parse("deducereturn")]      := deducereturn;
  TrussMe:-FEM:-m_CodegenOptions[parse("coercetypes")]       := coercetypes;
  TrussMe:-FEM:-m_CodegenOptions[parse("deducetypes")]       := deducetypes;
  TrussMe:-FEM:-m_CodegenOptions[parse("reduceanalysis")]    := reduceanalysis;
  TrussMe:-FEM:-m_CodegenOptions[parse("defaulttype")]       := defaulttype;
  TrussMe:-FEM:-m_CodegenOptions[parse("functionprecision")] := functionprecision;
  return NULL;
end proc: # SetCodegenOptions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetCodegenOptions := proc(
  field::string := "all",
  $)::anything;

  description "Get options for code generation.";

  if (field = "all") then
    return TrussMe:-FEM:-m_CodegenOptions;
  else
    return TrussMe:-FEM:-m_CodegenOptions[parse(field)];
  end if;
end proc: # GetCodegenOptions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export TranslateToMatlab := proc(
  expr_list::list,
  $)

  description "Convert a list of expressions <expr_list> into Matlab code.";

  local i, lang_fncs;

  # Define language functions
  lang_fncs := ["floor", "ceil", "round", "trunc", "erf"];

  # Define new language
  CodeGeneration:-LanguageDefinition:-Define(
    "NewMatlab", extend = "Matlab",
    seq(AddFunction(
      lang_fncs[i], anything::numeric, lang_fncs[i]
    ), i = 1..nops(lang_fncs))
  );

  # Generate code
  return CodeGeneration:-Translate(
    expr_list,
    language = "NewMatlab",
    op(op(eval(m_CodegenOptions))),
    output = string
  );
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ApplyIndent := proc(
  ind::string,
  str::string,
  $)::string;

  description "Apply indentation <ind> to string <str>.";

  local out;

  # Apply initial indentation
  out := cat(ind, str);

  # Apply indentation to each line
  out := StringTools:-SubstituteAll(out, "\n", cat("\n", ind));

  # Remove indentation from last line
  out := StringTools:-Take(out, length(out) - length(ind));

  # Remove indentation empty lines or lines with only spaces
  out := StringTools:-SubstituteAll(out, "\n" || ind || "\n", "\n\n");

  return out
end proc: # ApplyIndent

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateProperties := proc(
  data::list(symbol),
  {
  indent::string := "  "
  }, $)::string;

  description "Generate properties code from a list of data <data> and "
    "optional indentation <indent>.";

  local i, out;

  if (nops(data) > 0) then
    out := "";
    for i from 1 to nops(data) do
      convert(data[i], string);
      out := cat(out, indent, %, " = this.m_data.", %, ";\n");
    end do;
  else
    out := cat(indent, "% No properties\n");
  end if;
  return out;
end proc: # GenerateProperties

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateInputs := proc(
  vars::list(list({symbol, string})),
  {
  skipnull::boolean := true,
  indent::string    := "  "
  }, $)::string;

  description "Generate inputs code from a list of variables <vars> and "
    "optional indentation <indent> and skip null inputs <skipnull>.";

  local i, j, out;

  out := "";
  for j from 1 to nops(vars) do
    for i from 1 to nops(vars[j]) do
      if not (skipnull and type(vars[j][i], string)) then
        convert(vars[j][i], string);
        out := cat(out, indent, %, " = in_", j, "(", i, ");\n");
      end if;
    end do;
  end do;

  if evalb(out <> "") then
    return out;
  else
    return cat(indent, "% No inputs\n");
  end if;
end proc: # GenerateInputs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ExtractElements := proc(
  name::string,
  func::{list, Vector, Matrix, Array},
  dims::list(nonnegint),
  {
  skipnull::boolean := true,
  label::string     := "out",
  indent::string    := "  "
  }, $)::list, string;

  description "Extract elements for a n-dimensional function <func> with "
    "name <name>, dimensions <dims>, and optional veiling label <label>, "
    "indentation <indent> and skip null inputs <skipnull>.";

  local i, j, idx, lst, out, cur, tmp, str1, str2;

  lst := [];
  out := "";
  if has(map(type, func, numeric), false) then

    # Symbolic entries
    for i from 1 to mul(dims) do
      cur := [];
      idx := i - 1;
      for j from 1 to nops(dims) do
        cur := [op(cur), irem(idx, dims[j], 'idx') + 1];
      end do;
      tmp := func[op(cur)];
      if skipnull and evalb(tmp <> 0) then
        map(x -> (x, "_"), cur); str1 := cat(label, "_", op(%))[1..-2];
        lst := [op(lst), convert(str1, symbol) = tmp];
        map(x -> (x, ", "), cur); str2 := cat("", op(%))[1..-3];
        out := cat(out,
          indent, "out_", name, "(", str2, ") = ", str1, ";\n"
        );
      end if;
    end do;

  else

    # Numeric entries
    for i from 1 to mul(dims) do
      cur := [];
      idx := i - 1;
      for j from 1 to nops(dims) do
        cur := [op(cur), irem(idx, dims[j], 'idx') + 1];
      end do;
      tmp := func[op(cur)];
      if skipnull and evalb(tmp <> 0) then
        map(x -> (x, ", "), cur); str2 := cat("", op(%))[1..-3];
        out := cat(out,
          indent, "out_", name, "(", str2, ") = ", func[op(cur)], ";\n"
        );
      end if;
    end do;

  end if;
  return lst, out;
end proc: # ExtractElements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateHeader := proc(
  name::string,
  vars::list(list(symbol)),
  {
  skipthis::boolean := false,
  info::string      := "No description provided.",
  indent::string    := "  "
  }, $)::string;

  description "Generate a function header for a function <name> with "
    "variables <vars>, optional description <info> and indentation "
    "<indent> and skip class object input <skipthis>.";

  local i, out;

  out := cat("function out_", name, " = ", name, "( ");
  if skipthis and evalb(add(map(nops, vars)) = 0) then
    out := cat(out, "~");
  else
    out := cat(out, "this");
  end if;

  for i from 1 to nops(vars) do
    if evalb(nops(vars[i]) > 0) then
      out := cat(out, ", ", "in_", i);
    else
      out := cat(out, ", ", "~");
    end if;
  end do;
  return cat(out, " )\n", indent, "% ", info, "\n\n");
end proc: # GenerateHeader

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateElements := proc(
  func::{list, Vector, Matrix, Array},
  {
  indent::string := "  "
  }, $)::string;

  description "Generate code for elements <func> with optional indentation "
    "<indent>.";

  TrussMe:-FEM:-TranslateToMatlab(func);
  return TrussMe:-FEM:-ApplyIndent(indent,
    `if`(% = "", "% No body\n", cat(%%, %))
  );
end proc: # GenerateElements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateBody := proc(
  name::string,
  dims::list(nonnegint),
  {
  header::string     := "No header",
  properties::string := "No properties",
  inputs::string     := "No inputs",
  veils::string      := "No veils",
  elements::string   := "No elements",
  outputs::string    := "No outputs",
  indent::string     := "  ",
  typestr::string    := "zeros"
  }, $)::string;

  description "Generate code for function body for a function <name> with "
    "dimensions <dims>, optional header <header>, properties <properties>, "
    "inputs <inputs>, veils <veils>, elements <elements>, indentation "
    "<indent>, outputs <outputs> and vector type <typestr>.";

  local tmp;

  tmp := cat(seq(cat(dims[i], ", "), i = 1..nops(dims)-1), dims[-1]);
  if (nops(dims) <= 1) then
    tmp := cat(tmp, ", 1");
  end if;

  return cat(
    header,
    indent, "% Extract properties\n", properties, "\n",
    indent, "% Extract inputs\n",     inputs,     "\n",
    indent, "% Evaluate function\n",  elements,   "\n",
    indent, "% Store outputs\n",
    indent, "out_", name, " = ", typestr, "(", tmp, ");\n",
    outputs,
    "end % ", name, "\n"
  );
end proc: # GenerateBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export VectorToMatlab := proc(
  name::string,
  vars::list(list(symbol)),
  vec::Vector,
  {
  skipnull::boolean  := true,
  data::list(symbol) := [],
  info::string       := "No info",
  label::string      := "out",
  indent::string     := "  "
  }, $)::string;

  description "Translate the vector <vec> with variables <vars> into a "
    "Matlab function named <name> and return it as a string. The optional "
    "arguments and class properties <data>, function description <info>, "
    "veiling label <label>, and indentation string <indent>.";

  local header, properties, inputs, veils, elements, outputs, dims, lst,
    tmp_data, vec_inds, tmp_vars, tmp_vars_rm, tmp_vars_nl, typestr;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Generating vector function '%s'... ", name);
  end if;

  # Extract the function properties
  tmp_data := convert(data, set) intersect indets(vec, symbol);
  tmp_data := remove[flatten](j -> not evalb(j in tmp_data), data);
  properties := TrussMe:-FEM:-GenerateProperties(
    tmp_data, parse("indent") = indent
  );

  # Extract the function elements
  dims := [LinearAlgebra:-Dimension(vec)];
  lst, outputs := TrussMe:-FEM:-ExtractElements(
    name, vec, dims, parse("skipnull") = skipnull, parse("indent") = indent,
    parse("label") = label
    );

  # Extract the function inputs
  vec_inds    := indets(vec, symbol);
  tmp_vars    := map(i -> i intersect vec_inds, map(convert, vars, set));
  tmp_vars_rm := zip((i, j) -> remove[flatten](k -> not evalb(k in j), i), vars, tmp_vars);
  tmp_vars_nl := zip((i, j) -> map(k -> `if`(not evalb(k in j), "", k), i), vars, tmp_vars);
  inputs := TrussMe:-FEM:-GenerateInputs(
    tmp_vars_nl, parse("indent") = indent, parse("skipnull") = true
  );

  # Generate the method header
  header := TrussMe:-FEM:-GenerateHeader(
    name, tmp_vars_rm, parse("info") = info, parse("indent") = indent,
    parse("skipthis") = evalb(nops(tmp_data) = 0)
  );

  # Generate the elements
  elements := TrussMe:-FEM:-GenerateElements(lst);

  # Check if the vector is sparse
  if evalb(rtable_scanblock(vec, [], ':-NonZeros') < 0.5*dims[1]) and
     evalb(dims[1] > 5) then
    typestr := "sparse";
  else
    typestr := "zeros";
  end if;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  # Generate the generated code
  return TrussMe:-FEM:-GenerateBody(
    name, dims,
    parse("header")  = header, parse("properties") = properties,
    parse("inputs")  = inputs, parse("elements")   = elements,
    parse("indent")  = indent, parse("outputs")    = outputs,
    parse("typestr") = typestr
  );
end proc: # VectorToMatlab

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MatrixToMatlab := proc(
  name::string,
  vars::list(list(symbol)),
  mat::Matrix,
  {
  skipnull::boolean  := true,
  data::list(symbol) := [],
  info::string       := "No info",
  label::string      := "out",
  indent::string     := "  "
  }, $)::string;

  description "Translate the vector <vec> with variables <vars> into a "
    "Matlab function named <name> and return it as a string. The optional "
    "arguments and class properties <data>, function description <info>, "
    "veiling label <label>, and indentation string <indent>.";

  local header, properties, inputs, veils, elements, outputs, dims, lst,
    tmp_data, mat_inds, tmp_vars, tmp_vars_rm, tmp_vars_nl, typestr;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("Generating matrix function '%s'... ", name);
  end if;

  # Extract the function properties
  tmp_data := convert(data, set) intersect indets(mat, symbol);
  tmp_data := remove[flatten](j -> not evalb(j in tmp_data), data);
  properties := TrussMe:-FEM:-GenerateProperties(
    tmp_data, parse("indent") = indent
  );

  # Extract the function elements
  dims := [LinearAlgebra:-Dimensions(mat)];
  lst, outputs := TrussMe:-FEM:-ExtractElements(
    name, mat, dims, parse("skipnull") = skipnull, parse("indent") = indent,
    parse("label") = label
  );

  # Extract the function inputs
  mat_inds    := indets(mat, symbol);
  tmp_vars    := map(i -> i intersect mat_inds, map(convert, vars, set));
  tmp_vars_rm := zip((i, j) -> remove[flatten](k -> not evalb(k in j), i), vars, tmp_vars);
  tmp_vars_nl := zip((i, j) -> map(k -> `if`(not evalb(k in j), "", k), i), vars, tmp_vars);
  inputs := TrussMe:-FEM:-GenerateInputs(
    tmp_vars_nl, parse("indent") = indent, parse("skipnull") = true
  );

  # Generate the method header
  header := TrussMe:-FEM:-GenerateHeader(
    name, tmp_vars_rm, parse("info") = info, parse("indent") = indent,
    parse("skipthis") = evalb(nops(tmp_data) = 0)
  );

  # Generate the elements
  elements := TrussMe:-FEM:-GenerateElements(lst);

  # Check if the vector is sparse
  if evalb(rtable_scanblock(mat, [], ':-NonZeros') < 0.5*dims[1]*dims[2]) and
     evalb(dims[1]*dims[2] > 25) then
    typestr := "sparse";
  else
    typestr := "zeros";
  end if;

  if TrussMe:-FEM:-m_VerboseMode then
    printf("DONE\n");
  end if;

  # Store the results
  return TrussMe:-FEM:-GenerateBody(
    name, dims,
    parse("header")  = header, parse("properties") = properties,
    parse("inputs")  = inputs, parse("elements")   = elements,
    parse("indent")  = indent, parse("outputs")    = outputs,
    parse("typestr") = typestr
  );
end proc: # MatrixToMatlab

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetVeilSubs := proc(
  veil::{Vector(algebraic), list(algebraic)},
  $)

  description "Get the substitutions to transform veils <veil> from "
    "(v[j])(f) to v_j.";

  local v, rm_v_deps, i, v_tmp;

  # Store veiling variables
  if not type(veil, list) then
    v := convert(veil, list);
  else
    v := veil;
  end if;

  # Compute transformations
  if (nops(v) > 0) then
    # Transform veil variables (v[i])(f) -> v_i
    v_tmp := Array(v);
    #for i from 1 to nops(v_tmp) do
    for i from 1 to op(2, ArrayDims(v_tmp)) do
      # (v[i])(f) -> v[i]
      if type(v_tmp[i], function) then
        v_tmp[i] := op(0, v_tmp[i]);
      end if;
      # v[i] -> v_i
      if type(v_tmp[i], indexed) then
        v_tmp[i] := convert(
          cat(op(0, v_tmp[i]), "_", op(1..-1, v_tmp[i])), symbol
        );
      end if;
    end do;
    rm_v_deps := v =~ convert(v_tmp, list);
  else
    rm_v_deps := [];
  end if;

  # Return output
  return rm_v_deps;
end proc: # GetVeilSubs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateConstructor := proc(
  name::string,
  {
  data::list(symbol = algebraic) := [],
  info::string                   := "Class constructor.",
  indent::string                 := "  "
  }, $)::string;

  description "Generate a constructor for a system named <name> system data "
    "<data>, description <info> and indentation <indent>.";

  return cat("function this = ", name, "( varargin )\n",
    TrussMe:-FEM:-ApplyIndent(
      indent, cat(
      "% ", info, "\n",
      "\n",
      "% User data\n",
      # No input arguments
      "if (nargin == 0)\n",
      `if`(evalb(nops(data) > 0),
      cat~(op(cat~(
        indent, "data.", convert~(lhs~(data), string),
        " = ", convert~(rhs~(data), string), ";\n"
      ))),
      cat(indent, "data = [];\n")
      ),
      # Struct input argument
      "elseif (nargin == 1 && isstruct(varargin{1}))\n",
      cat~(op(cat~(
        indent, "data = varargin{1};\n"
      ))),
      # Many input arguments
      `if`(evalb(nops(data) > 0),
      cat("elseif (nargin == ", nops(data), ")\n",
      cat~(op(cat~(
        indent, "data.", convert~(lhs~(data), string),
        " = varargin{", convert~([seq(1..nops(data))], string), "};\n"
      )))), ""),
      # Wrong number of input arguments
      "else\n",
      indent, "error('wrong number of input arguments.');\n",
      "end\n"
      "\n",
      "% Call superclass constructor\n",
      "this = this@TrussMe.System(data);\n"
    )),
    "end % ", name, "\n");
end proc: # GenerateConstructor

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SystemToMatlab := proc(
  name::string,
  fem::FEM,
  {
  vars::list(symbol)             := [],
  data::list(symbol = algebraic) := [],
  info::string                   := "No class description provided.",
  label::string                  := "out",
  indent::string                 := "  "
  }, $)::string;

  description "Generate a FEM system <fem> with name <name>, with optional "
    "data <data>, description <info>, output label <label> and indentation "
    "<indent>.";

  local i, bar, rm_v_deps, data_vars, x, v, props_str, check_str;

  # Function utilities strings
  i   := indent;
  bar := "% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";

  # Get veiling variables substitution
  rm_v_deps := TrussMe:-FEM:-GetVeilSubs(lhs~(fem["veils"]));

  # Retrieve variables and data
  x := vars;
  v := subs(op(rm_v_deps), convert(lhs~(fem["veils"]), list));
  if evalb(nops(data) > 0) then
    data_vars := lhs~(data);
  else
    data_vars := [];
  end if;

  # Prepare a string to check to detect if the system has been solved
  if fem["solved"] then
    props_str := "";
    check_str := "";
  else
    props_str := cat(
      i, "properties (SetAccess = protected, Hidden = true)\n",
      i, i, "m_solved = ", fem["solved"], ";\n",
      i, "end\n",
      i, "%\n"
    );
    check_str := cat(
      i, "\n\n",
      i, "% Check if the system has been solved\n",
      i, "if ~this.m_solved\n",
      i, i, "error('system has not been solved.');\n",
      i, "end"
    );
  end if;

  # Return output string
  return cat(
    "% +--------------------------------------------------------------------------+\n",
    "% | 'TrussMe' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023 |\n",
    "% | Current version authors:                                                 |\n",
    "% |   Davide Stocco and Matteo Larcher.                                      |\n",
    "% +--------------------------------------------------------------------------+\n",
    "\n",
    "% Matlab generated code for system: ", name, "\n",
    "% This file has been automatically generated by TrussMe.\n",
    "% DISCLAIMER: If you need to edit it, do it wisely!\n",
    "\n",
    "classdef ", name, " < TrussMe.System\n",
    i, "%\n",
    i, "% ", info, "\n",
    i, "%\n",
    props_str,
    i, "methods\n",
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      GenerateConstructor(
        name,
        parse("data")   = data,
        parse("info")   = "Class constructor.",
        parse("indent") = i
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-MatrixToMatlab(
        "K", [x, []], subs(op(rm_v_deps), fem["K"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the stiffness matrix K."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-MatrixToMatlab(
        "K_ff", [x, []], subs(op(rm_v_deps), fem["K_ff"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the stiffness matrix K_ff."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-MatrixToMatlab(
        "K_fs", [x, []], subs(op(rm_v_deps), fem["K_fs"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the stiffness matrix K_fs."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-MatrixToMatlab(
        "K_sf", [x, []], subs(op(rm_v_deps), fem["K_sf"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the stiffness matrix K_sf."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-MatrixToMatlab(
        "K_ss", [x, []], subs(op(rm_v_deps), fem["K_ss"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the stiffness matrix K_ss."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "d", [x, v], subs(op(rm_v_deps), fem["d"]),
        parse("data") = data_vars,
        parse("info") = cat("Evaluate the deformation vector d.", check_str)
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "d_f", [x, v], subs(op(rm_v_deps), fem["d_f"]),
        parse("data") = data_vars,
        parse("info") = cat("Evaluate the deformation vector d_f.", check_str)
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "d_s", [x, []], subs(op(rm_v_deps), fem["d_s"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the deformation vector d_s."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "f", [x, v], subs(op(rm_v_deps), fem["f"]),
        parse("data") = data_vars,
        parse("info") = cat("Evaluate the force vector f.", check_str)
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "f_f", [x, []], subs(op(rm_v_deps), fem["f_f"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the force vector f_f."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "f_s", [x, v], subs(op(rm_v_deps), fem["f_s"]),
        parse("data") = data_vars,
        parse("info") = cat("Evaluate the force vector f_s.", check_str)
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "f_r", [x, []], subs(op(rm_v_deps), fem["f_r"]),
        parse("data") = data_vars,
        parse("info") = "Evaluate the force vector f_r."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "perm", [], subs(op(rm_v_deps), convert(fem["perm"], Vector)),
        parse("data") = data_vars,
        parse("info") = "Evaluate the permutation vector."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "unperm", [], subs(op(rm_v_deps), convert(fem["unperm"], Vector)),
        parse("data") = data_vars,
        parse("info") = "Evaluate the unpermutation vector."
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    TrussMe:-FEM:-ApplyIndent(
      cat(i, i),
      TrussMe:-FEM:-VectorToMatlab(
        "v", [x], subs(op(rm_v_deps), convert(rhs~(fem["veils"]), Vector)),
        parse("data")  = data_vars,
        parse("label") = fem["label"],
        parse("info")  = cat("Evaluate the veiling vector.", check_str)
    )),
    i, i, "%\n",
    i, i, bar,
    i, i, "%\n",
    i, "end\n",
    "end % ", name, "\n",
    "\n",
    "% That's All Folks!\n"
  );
end proc: # SystemToMatlab

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GenerateMatlabCode := proc(
  name::string,
  fem::FEM,
  {
  path::string                   := "./",
  vars::list(symbol)             := [],
  data::list(symbol = algebraic) := [],
  info::string                   := "No class description provided.",
  label::string                  := "out",
  indent::string                 := "  "
  }, $)

  description "Generate Matlab code for the FEM system <fem> with name <name>, "
    "with optional data <data>, description <info>, output label <label> and "
    "indentation <indent>, variables <vars> and output path <path>.";

  TrussMe:-FEM:-GenerateFile(
    cat(path, name, ".m"),
    TrussMe:-FEM:-SystemToMatlab(
      name, fem,
      parse("vars")   = vars,
      parse("data")   = data,
      parse("info")   = info,
      parse("label")  = label,
      parse("indent") = indent
  ));
  return NULL;
end proc: # GenerateMatlabCode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
