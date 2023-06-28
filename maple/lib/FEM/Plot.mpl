# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                             ____  _       _                                 #
#                            |  _ \| | ___ | |_                               #
#                            | |_) | |/ _ \| __|                              #
#                            |  __/| | (_) | |_                               #
#                            |_|   |_|\___/ \__|                              #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco  (University of Trento)
#   Matteo Larcher (University of Trento)
#
# License: BSD 3-Clause License

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ObjectColor := proc(
  obj::{NODE, SUPPORT, ELEMENT, LOAD},
  $)::string;

  description "Return the color of the object <obj>.";

  if type(obj, SUPPORT) then
    return TrussMe:-FEM:-m_SupportColor;
  elif type(obj, NODE) and not type(obj, SUPPORT) then
    return TrussMe:-FEM:-m_NodeColor;
  elif type(obj, ELEMENT) and evalb(nops(obj["nodes"]) = 2) then
    return TrussMe:-FEM:-m_ElementColor;
  elif type(obj, ELEMENT) and evalb(nops(obj["nodes"]) > 2) then
    return TrussMe:-FEM:-m_ShellColor;
  elif type(obj, LOAD) then
    return TrussMe:-FEM:-m_ForceColor;
  else
    error("<obj> must be a node, support, element or load.");
  end if;
end proc: # ObjectColor

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export StructureGraph := proc(
  fem::FEM,
  {
  id::boolean   := false,
  disp::boolean := true
  }, $)::function;

  description "Return the graph of the structure given a list of nodes "
    "<nodes> and elements <elements>.";

  local nodes, elements, loads, obj, i, vertex_name, vertex_id, vertex_color, G;

  if m_VerboseMode then
    printf("Checking structure connections... ");
  end if;

  # Extract nodes, elements and loads from <fem> object
  nodes    := fem["nodes"];
  elements := fem["elements"];
  loads    := fem["loads"];

  # Retrieve objects names, ids and colors
  vertex_name := [
    seq(obj["name"], obj = nodes),
    seq(obj["name"], obj = elements),
    seq(obj["name"], obj = loads)
  ];
  vertex_id := [
    seq(obj["id"], obj = nodes),
    seq(obj["id"], obj = elements),
    seq(obj["id"], obj = loads)
  ];
  vertex_color := [
    seq(TrussMe:-FEM:-ObjectColor(obj), obj = nodes),
    seq(TrussMe:-FEM:-ObjectColor(obj), obj = elements),
    seq(TrussMe:-FEM:-ObjectColor(obj), obj = loads)
  ];

  # Build graph
  G := GraphTheory:-Graph(vertex_id);
  GraphTheory:-HighlightVertex(G, vertex_id, vertex_color);

  # Connect nodes with elements
  for obj in elements do
    for i from 1 to nops(obj["nodes"]) do
      GraphTheory:-AddEdge(G, {obj["id"], obj["nodes"][i]});
    end do;
  end do;

  # Connect nodes with loads
  for obj in loads do
    GraphTheory:-AddEdge(G, {obj["id"], obj["node"]});
  end do;

  # Relabel with vertex names
  if not id then
    G := GraphTheory:-RelabelVertices(G, vertex_name);
  end if;

  if TrussMe:-FEM:-m_VerboseMode then
    printf(" DONE\n");
  end if;

  # Plot graph
  if disp then
    print(plots:-display(GraphTheory:-DrawGraph(G, layout = tree)));
  end if;

  # Return graph
  return G;
end proc: # StructureGraph

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DrawFrame := proc(
  frame::FRAME,
  {
  data::{list(`=`), set(`=`)} := [],
  scaling::numeric            := 1.0,
  colors::list(string)        := ["Red", "Green", "Blue"]
  }, $)::function;

  description "Draw a reference frame <frame> given a list of substitution "
  "data <data>, an axes scaling factor <scaling>, and axes colors <colors>.";

  local tmp, p, x, y, z:

  if evalb(nops(colors) <> 3) then
    error("<colors> must be a list of 3 elements.");
  end if;

  tmp := subs(op(data), frame);
  p := [TrussMe:-FEM:-CompXYZ(TrussMe:-FEM:-Origin(tmp))]:
  x := [TrussMe:-FEM:-CompXYZ(TrussMe:-FEM:-UvecX(tmp))]:
  y := [TrussMe:-FEM:-CompXYZ(TrussMe:-FEM:-UvecY(tmp))]:
  z := [TrussMe:-FEM:-CompXYZ(TrussMe:-FEM:-UvecZ(tmp))]:

  return plots:-display(
    plots:-arrow(p, scaling*x, parse("shape") = cylindrical_arrow, parse("color") = colors[1]),
    plots:-arrow(p, scaling*y, parse("shape") = cylindrical_arrow, parse("color") = colors[2]),
    plots:-arrow(p, scaling*z, parse("shape") = cylindrical_arrow, parse("color") = colors[3])
  );
end proc: # DrawFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotNode := proc(
  p::{list(algebraic), Vector(algebraic)},
  {
  data::{list(`=`), set(`=`)} := [],
  token::symbol               := TrussMe:-FEM:-m_NodeToken,
  color::string               := TrussMe:-FEM:-m_NodeColor
  }, $)::function;

  description "Plot the node (or support) at point <p> given a list or set of "
    "data for substitution <data>, a display token <token> and a display color "
    "<color>.";

  local p_tmp;

  if type(p, list) and evalb(nops(p) = 3) then
    p_tmp := subs(op(data), p);
  elif type(p, Vector) and evalb(LinearAlgebra:-Dimension(p) = 4) then
    p_tmp := convert(subs(op(data), p), list)[1..3];
  else
    error("<p> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  return plots:-display(
    plottools:-point(p_tmp, parse("symbol") = token, parse("symbolsize") = 15),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotElement := proc(
  p_1::{list(algebraic), Vector(algebraic)},
  p_2::{list(algebraic), Vector(algebraic)},
  {
  data::{list(`=`), set(`=`)} := [],
  color::string               := TrussMe:-FEM:-m_ElementColor
  }, $)::function;

  description "Plot the element from point <p_1> and <p_2> given a list or set "
    "of data for substitution <data> and a display color <color>.";

  local p_1_tmp, p_2_tmp;

  if type(p_1, list) and evalb(nops(p_1) = 3) then
    p_1_tmp := subs(op(data), p_1);
  elif type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 4) then
    p_1_tmp := convert(subs(op(data), p_1), list)[1..3];
  else
    error("<p_1> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp := subs(op(data), p_2);
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := convert(subs(op(data), p_2), list)[1..3];
  else
    error("<p_2> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  return plots:-display(
    plottools:-line(p_1_tmp, p_2_tmp, parse("thickness") = 6),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedElement := proc(
  p_1::{list(algebraic), Vector(algebraic)},
  p_2::{list(algebraic), Vector(algebraic)},
  d_1::{list(algebraic), Vector(algebraic)},
  d_2::{list(algebraic), Vector(algebraic)},
  {
  data::{list(`=`), set(`=`)} := [],
  scaling::nonnegative        := 1.0,
  color::string               := TrussMe:-FEM:-m_ElementColor
  }, $)::function;

  description "Plot the element from diplacements <d_1> and <d_2> given a list "
    "or set of data for substitution <data> and a display color <color>.";

  local d_1_tmp, d_2_tmp, p_1_tmp, p_2_tmp, p1p2, p1p2_unit, x, y, z,
    deformed_element;

  if type(d_1, list) and evalb(nops(d_1) = 6) then
    d_1_tmp := convert(subs(op(data), d_1), Vector);
  elif type(d_1, Vector) and evalb(LinearAlgebra:-Dimension(d_1) = 6) then
    d_1_tmp := subs(op(data), d_1);;
  else
    error("<d_1> must be a list or a vector of 6 elements.");
  end if;

  if type(d_2, list) and evalb(nops(d_2) = 6) then
    d_2_tmp := convert(subs(op(data), d_2), Vector);
  elif type(d_2, Vector) and evalb(LinearAlgebra:-Dimension(d_2) = 6) then
    d_2_tmp := subs(op(data), d_2);
  else
    error("<d_2> must be a list or a vector of 6 elements.");
  end if;

  if type(p_1, list) and evalb(nops(p_1) = 3) then
    p_1_tmp := <convert(subs(op(data), p_1), Vector), 1>;
  elif type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 4) then
    p_1_tmp := subs(op(data), p_1);
  else
    error("<p_1> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp := <convert(subs(op(data), p_2), Vector), 1>;
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := subs(op(data), p_2);
  else
    error("<p_2> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  p1p2 := p_2_tmp - p_1_tmp;
  p1p2_unit := try abs(p1p2 /~ TrussMe:-FEM:-Norm2(p1p2[1..3])) catch: [1,0,0] end try;

  # Calculate intermediate points through interpolation
  x := (d_1_tmp[1]*(1-xi) + d_2_tmp[1]*xi) * p1p2_unit[1] +
       (d_1_tmp[1]*(1-3*xi^2+2*xi^3) + d_2_tmp[1]*(3*xi^2-2*xi^3)) *
       TrussMe:-FEM:-Norm2(p1p2_unit[2..3]) +
       (d_1_tmp[6]*(xi-xi^2)*(1-xi) + d_2_tmp[6]*(xi^2-xi^3)*xi) *
       TrussMe:-FEM:-Norm2(p1p2_unit[2..3]) -
       (d_1_tmp[5]*(xi-xi^2)*(1-xi) + d_2_tmp[5]*(xi^2-xi^3)*xi) *
       TrussMe:-FEM:-Norm2(p1p2_unit[2..3]);
  y := (d_1_tmp[2]*(1-xi) + d_2_tmp[2]*xi) * p1p2_unit[2] +
       (d_1_tmp[2]*(1-3*xi^2+2*xi^3) + d_2_tmp[2]*(3*xi^2-2*xi^3)) *
       TrussMe:-FEM:-Norm2(p1p2_unit[[1,3]]) -
       (d_1_tmp[6]*(xi-xi^2)*(1-xi) + d_2_tmp[6]*(xi^2-xi^3)*xi) *
       TrussMe:-FEM:-Norm2(p1p2_unit[[1,3]]) +
       (d_1_tmp[4]*(xi-xi^2)*(1-xi) + d_2_tmp[4]*(xi^2-xi^3)*xi) *
       TrussMe:-FEM:-Norm2(p1p2_unit[[1,3]]);
  z := (d_1_tmp[3]*(1-xi) + d_2_tmp[3]*xi) * p1p2_unit[3] +
       (d_1_tmp[3]*(1-3*xi^2+2*xi^3) + d_2_tmp[3]*(3*xi^2-2*xi^3)) *
       TrussMe:-FEM:-Norm2(p1p2_unit[1..2]) -
       (d_1_tmp[4]*(xi-xi^2)*(1-xi) + d_2_tmp[4]*(xi^2-xi^3)*xi) *
       TrussMe:-FEM:-Norm2(p1p2_unit[1..2]) +
       (d_1_tmp[5]*(xi-xi^2)*(1-xi) + d_2_tmp[5]*(xi^2-xi^3)*(xi)) *
       TrussMe:-FEM:-Norm2(p1p2_unit[1..2]);

  deformed_element := scaling *~ <x, y, z, 0> + p_1_tmp + p1p2 *~ xi;

  return plots:-display(
    plots:-spacecurve(
      convert(deformed_element, list)[1..3], xi = 0..1, parse("thickness") = 6
    ),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotDeformedElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotLoad := proc(
  p_1::{list(algebraic), Vector(algebraic)},
  p_2::{list(algebraic), Vector(algebraic)},
  {
  data::{list(`=`), set(`=`)} := [],
  color::string               := TrussMe:-FEM:-m_LoadColor,
  scaling::nonnegative        := 1.0
  }, $)::function;

  description "Plot the load arrow from point <p_1> and <p_2> given a list or "
    "set of data for substitution <data> and a display color <color>.";

  local p_1_tmp, p_2_tmp, wb, wh, hh;

  if type(p_1, list) and evalb(nops(p_1) = 3) then
    p_1_tmp :=  subs(op(data), p_1);
  elif type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 4) then
    p_1_tmp := convert( subs(op(data), p_1), list)[1..3];
  else
    error("<p_1> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp :=  subs(op(data), p_2);
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := convert( subs(op(data), p_2), list)[1..3];
  else
    error("<p_2> must be a list of 3 elements or a vector of 4 elements.");
  end if;

  wb := max(0.02, 0.05 * scaling);
  wh := max(0.04, 0.10 * scaling);
  hh := min(0.1, try
      (0.1 * scaling)/TrussMe:-FEM:-Norm2(p_2_tmp - p_1_tmp) catch: infinity end try
  );
  return plots:-display(
    plottools:-arrow(p_1_tmp, p_2_tmp, wb, wh, hh, cylindrical_arrow),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotLoad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotStructure := proc(
  fem::FEM,
  {
  data::{list(`=`), set(`=`)}             := [],
  frame_scaling::{numeric, list(numeric)} := 0.0,
  load_scaling::numeric                   := 1.0
  }, $)::{function, list(function)};

  description"Plot the undeformed <fem> structure given a list or set of "
    "substitution data <data>, a frame scaling factor <frame_scaling> or "
    "<nodes_scaling, elements_scaling>, and a loads scaling factor <scaling>.";

  local nodes, elements, loads, tmp, tmp_i, n, i, j, k, p_1, p_2,
    nodes_scaling, elements_scaling, disp_nodes_frames, disp_elements_frames,
    disp_nodes, disp_elements, disp_loads;

  # Extract the nodes, elements and loads from <fem> object
  nodes    := fem["nodes"];
  elements := fem["elements"];
  loads    := fem["loads"];

  # Retrieve frame scaling factors
  if type(frame_scaling, numeric) then
    nodes_scaling := frame_scaling;
    elements_scaling := frame_scaling;
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 1) then
    nodes_scaling := frame_scaling[1];
    elements_scaling := frame_scaling[1];
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 2) then
    nodes_scaling := frame_scaling[1];
    elements_scaling := frame_scaling[2];
  else
    error("<frame_scaling> or <nodes_scaling, elements_scaling> must be a list "
      "of 1 or 2 elements.");
  end if;

  # Draw nodes frames
  if evalb(nodes_scaling <> 0) then
    disp_nodes_frames := [seq(i, i = 1..nops(nodes))];
    for i from 1 to nops(nodes) do
      disp_nodes_frames[i] := TrussMe:-FEM:-DrawFrame(
        nodes[i]["frame"], parse("scaling") = nodes_scaling
      );
    end do;
  else
    disp_nodes_frames := [];
  end if;

  # Draw elements frames
  n := nops(elements);
  if evalb(elements_scaling <> 0) then
    disp_elements_frames := [seq(i, i = 1..n)];
    for i from 1 to n do
      disp_elements_frames[i] := TrussMe:-FEM:-DrawFrame(
        elements[i]["frame"], parse("scaling") = elements_scaling
      );
    end do;
  else
    disp_elements_frames := [];
  end if;

  # Plot nodes
  disp_nodes := [seq(i, i = 1..nops(nodes))];
  for i from 1 to nops(nodes) do
    p_1 := <nodes[i]["coordinates"], 1>;
    disp_nodes[i] := TrussMe:-FEM:-PlotNode(
      convert(p_1, Vector), parse("data")  = data,
      parse("color") = `if`(
        type(nodes[i], SUPPORT), TrussMe:-FEM:-m_SupportColor, TrussMe:-FEM:-m_NodeColor
      ),
      parse("token") = `if`(
        type(nodes[i], SUPPORT), TrussMe:-FEM:-m_SupportToken, TrussMe:-FEM:-m_NodeToken
      )
    );
  end do;

  # Plot elements
  disp_elements := [seq(i, i = 1..add(
    map(i -> nops(i["nodes"])*(nops(i["nodes"])-1)/2, elements)
  ))];
  tmp_i := 0;
  for i from 1 to n do
    tmp := nops(elements[i]["nodes"])-1;
    for j from 1 to tmp do
      for k from j to tmp do
        tmp_i := tmp_i + 1;
        TrussMe:-FEM:-GetObjById(nodes, elements[i]["nodes"][j], parse("position") = true);
        p_1 := <nodes[%]["coordinates"], 1>;
        TrussMe:-FEM:-GetObjById(nodes, elements[i]["nodes"][k+1], parse("position") = true);
        p_2 := <nodes[%]["coordinates"], 1>;
        disp_elements[tmp_i] := TrussMe:-FEM:-PlotElement(
          convert(p_1, Vector), convert(p_2, Vector),
          parse("data")  = data,
          parse("color") = `if`(evalb(nops(elements[i]["nodes"]) = 2),
            TrussMe:-FEM:-m_ElementColor, TrussMe:-FEM:-m_ShellColor)
        );
      end do;
    end do;
  end do;

  # Plot loads
  if evalb(load_scaling <> 0) then
    disp_loads := [seq(i, i = 1..2*nops(loads))];
    for i from 1 to nops(loads) do
      j := TrussMe:-FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
      p_2 := <nodes[j]["coordinates"], 1>;

      # Plot forces
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      end if;
      disp_loads[2*i-1] := TrussMe:-FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe:-FEM:-m_ForceColor
      );

      # Plot moments
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      end if;
      disp_loads[2*i] := TrussMe:-FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe:-FEM:-m_MomentColor
      );
    end do;
  else
    disp_loads := [];
  end if;

  # Plot structure
  return plots:-display(
    [op(disp_nodes_frames), op(disp_elements_frames),
     op(disp_nodes), op(disp_elements), op(disp_loads)],
    parse("axes")    = boxed,
    parse("scaling") = constrained,
    parse("labels")  = ['x', 'y', 'z']
  );
end proc: # PlotStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedStructure := proc(
  fem::FEM,
  {
  data::{list(`=`), set(`=`)}             := [],
  interpolate::boolean                    := true,
  frame_scaling::{numeric, list(numeric)} := 0.0,
  deformation_scaling::numeric            := 1.0,
  load_scaling::numeric                   := 1.0
  }, $)::{function, list(function)};

  description "Plot the deformed <fem> structure given a list or set of "
    "substitution data <data>, a frame scaling factor <frame_scaling> or "
    "<nodes_scaling, elements_scaling>, a loads scaling factor <load_scaling>, "
    "and a deformation magnification factor <deformation_scaling>.";

  local nodes, elements, loads, n, i, j, k, tmp, tmp_i, p_1, p_2, d_1, d_2,
    nodes_scaling, elements_scaling, disp_nodes_frames, disp_elements_frames,
    disp_nodes, disp_elements, disp_loads, data_tmp;

  # Extract the nodes, elements and loads from <fem> object
  nodes    := fem["nodes"];
  elements := fem["elements"];
  loads    := fem["loads"];

  # Add data evaluated and reversed veils to data substitution
  data_tmp := [op(subs(op(data), ListTools:-Reverse(fem["veils"]))), op(data)];

  # Retrieve frame scaling factors
  if type(frame_scaling, numeric) then
    nodes_scaling := frame_scaling;
    elements_scaling := frame_scaling;
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 1) then
    nodes_scaling := frame_scaling[1];
    elements_scaling := frame_scaling[1];
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 2) then
    nodes_scaling := frame_scaling[1];
    elements_scaling := frame_scaling[2];
  else
    error("<frame_scaling> or <nodes_scaling, elements_scaling> must be a list "
      "of 1 or 2 elements.");
  end if;

  # Draw nodes frames
  if evalb(nodes_scaling <> 0) then
    disp_nodes_frames := [seq(i, i = 1..nops(nodes))];
    for i from 1 to nops(nodes) do
      # TODO: translate and rotate according to node deformations
      disp_nodes_frames[i] := TrussMe:-FEM:-DrawFrame(
        nodes[i]["frame"], parse("scaling") = nodes_scaling
      );
    end do;
  else
    disp_nodes_frames := [];
  end if;

  # Draw elements frames
  n := nops(elements);
  if evalb(elements_scaling <> 0) then
    disp_elements_frames := [seq(i, i = 1..n)];
    for i from 1 to n do
      # TODO: translate and rotate according to node deformations
      disp_elements_frames[i] := TrussMe:-FEM:-DrawFrame(
        elements[i]["frame"], parse("scaling") = elements_scaling
      );
    end do;
  else
    disp_elements_frames := [];
  end if;

  # Plot deformed nodes
  disp_nodes := [seq(i, i = 1..nops(nodes))];
  for i from 1 to nops(nodes) do
    p_1 := convert(
      <nodes[i]["coordinates"] +
        deformation_scaling *~ nodes[i]["frame"][1..3, 1..3].nodes[i]["output_displacements"][1..3], 1>,
      Vector
    );
    disp_nodes[i] := TrussMe:-FEM:-PlotNode(
      p_1, parse("data")  = data_tmp,
      parse("color") = `if`(
        type(nodes[i], SUPPORT), TrussMe:-FEM:-m_SupportColor, TrussMe:-FEM:-m_NodeColor
      ),
      parse("token") = `if`(
        type(nodes[i], SUPPORT), TrussMe:-FEM:-m_SupportToken, TrussMe:-FEM:-m_NodeToken
      )
    );
  end do;

  # Plot deformed elements
  disp_elements := [seq(i, i = 1..add(
    map(i -> nops(i["nodes"])*(nops(i["nodes"])-1)/2, elements)
  ))];
  tmp_i := 0;
  for i from 1 to n do
    tmp := nops(elements[i]["nodes"]) - 1;
    for j from 1 to tmp do
      for k from j to tmp do
        tmp_i := tmp_i + 1;
        if not interpolate then
          TrussMe:-FEM:-GetObjById(nodes, elements[i]["nodes"][j], parse("position") = true);
          p_1 := <nodes[%]["coordinates"] +
            nodes[%]["frame"][1..3, 1..3].(deformation_scaling *~ nodes[%]["output_displacements"][1..3]), 1>;
          TrussMe:-FEM:-GetObjById(nodes, elements[i]["nodes"][k+1], parse("position") = true);
          p_2 := <nodes[%]["coordinates"] +
            nodes[%]["frame"][1..3, 1..3].(deformation_scaling *~ nodes[%]["output_displacements"][1..3]), 1>;
          disp_elements[tmp_i] := TrussMe:-FEM:-PlotElement(
            convert(p_1, Vector), convert(p_2, Vector), parse("data") = data,
            parse("color") = `if`(evalb(nops(elements[i]["nodes"]) = 2),
              TrussMe:-FEM:-m_ElementColor, TrussMe:-FEM:-m_ShellColor)
          );
        else
          TrussMe:-FEM:-GetObjById(nodes, elements[i]["nodes"][j], parse("position") = true);
          p_1 := <nodes[%]["coordinates"], 1>;
          d_1 := <nodes[%%]["frame"][1..3, 1..3].nodes[%%]["output_displacements"][1..3],
            nodes[%%]["frame"][1..3, 1..3].nodes[%%]["output_displacements"][4..6]>;
          TrussMe:-FEM:-GetObjById(nodes, elements[i]["nodes"][k+1], parse("position") = true);
          p_2 := <nodes[%]["coordinates"], 1>;
          d_2 := <nodes[%%]["frame"][1..3, 1..3].nodes[%%]["output_displacements"][1..3],
            nodes[%%]["frame"][1..3, 1..3].nodes[%%]["output_displacements"][4..6]>;
          disp_elements[tmp_i] := TrussMe:-FEM:-PlotDeformedElement(
            convert(p_1, Vector), convert(p_2, Vector),
            convert(d_1, Vector), convert(d_2, Vector),
            parse("data")    = data_tmp,
            parse("scaling") = deformation_scaling,
            parse("color")   = `if`(evalb(nops(elements[i]["nodes"]) = 2),
              TrussMe:-FEM:-m_ElementColor, TrussMe:-FEM:-m_ShellColor)
          );
        end if;
      end do;
    end do;
  end do;

  # Plot loads
  if evalb(load_scaling <> 0) then
    disp_loads := [seq(i, i = 1..2*nops(loads))];
    for i from 1 to nops(loads) do
      j := TrussMe:-FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
      p_2 := <nodes[j]["coordinates"] +
        deformation_scaling *~ nodes[j]["output_displacements"][1..3], 1>;

      # Plot forces
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      end if;
      disp_loads[2*i-1] := TrussMe:-FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data_tmp,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe:-FEM:-m_ForceColor
      );

      # Plot moments
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      end if;
      disp_loads[2*i] := TrussMe:-FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data_tmp,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe:-FEM:-m_MomentColor
      );
    end do;
  else
    disp_loads := [];
  end if;

  # Plot deformed structure
  return plots:-display(
    [op(disp_nodes_frames), op(disp_elements_frames),
     op(disp_nodes), op(disp_elements), op(disp_loads)],
    parse("axes")    = boxed,
    parse("scaling") = constrained,
    parse("labels")  = ['x', 'y', 'z']
  );
end proc: # PlotDeformedStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
