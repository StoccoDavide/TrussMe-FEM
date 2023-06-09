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
    return TrussMe_FEM:-m_SupportColor;
  elif type(obj, NODE) and not type(obj, SUPPORT) then
    return TrussMe_FEM:-m_NodeColor;
  elif type(obj, ELEMENT) then
    return TrussMe_FEM:-m_ElementColor;
  elif type(obj, LOAD) then
    return TrussMe_FEM:-m_ForceColor;
  else
    error "invalid object detected.";
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
    printf("Checking structure connections...");
  end if;

  # Extract the nodes, elements and loads from <fem> object
  nodes    := fem["nodes"];
  elements := fem["elements"];
  loads    := fem["loads"];

  # Retrieve the objects names, ids and colors
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
    seq(TrussMe_FEM:-ObjectColor(obj), obj = nodes),
    seq(TrussMe_FEM:-ObjectColor(obj), obj = elements),
    seq(TrussMe_FEM:-ObjectColor(obj), obj = loads)
  ];

  # Build the graph
  G := GraphTheory:-Graph(vertex_id);
  GraphTheory:-HighlightVertex(G, vertex_id, vertex_color);

  # Connect the nodes with the elements
  for obj in elements do
    GraphTheory:-AddEdge(G, {obj["id"], obj["node_1"]});
    GraphTheory:-AddEdge(G, {obj["id"], obj["node_2"]});
  end do;

  # Connect the nodes with the loads
  for obj in loads do
    GraphTheory:-AddEdge(G, {obj["id"], obj["node"]});
  end do;

  # Relabel with vertex names
  if not id then
    G := GraphTheory:-RelabelVertices(G, vertex_name);
  end if;

  if m_VerboseMode then
    printf("\tDONE\n");
  end if;

  # Plot the graph
  if disp then
    print(plots:-display(GraphTheory:-DrawGraph(G, layout = tree)));
  end if;

  # Return the graph
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
    error("<colors> expects a list of 3 elements.");
  end if;

  tmp := subs(op(data), frame);
  p := [TrussMe_FEM:-CompXYZ(TrussMe_FEM:-Origin(tmp))]:
  x := [TrussMe_FEM:-CompXYZ(TrussMe_FEM:-UvecX(tmp))]:
  y := [TrussMe_FEM:-CompXYZ(TrussMe_FEM:-UvecY(tmp))]:
  z := [TrussMe_FEM:-CompXYZ(TrussMe_FEM:-UvecZ(tmp))]:

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
  token::symbol               := TrussMe_FEM:-m_NodeToken,
  color::string               := TrussMe_FEM:-m_NodeColor
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
    error("invalid point vector detected.");
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
  color::string               := TrussMe_FEM:-m_ElementColor
  }, $)::function;

  description "Plot the element from point <p_1> and <p_2> given a list or set "
    "of data for substitution <data> and a display color <color>.";

  local p_1_tmp, p_2_tmp;

  if type(p_1, list) and evalb(nops(p_1) = 3) then
    p_1_tmp := subs(op(data), p_1);
  elif type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 4) then
    p_1_tmp := convert(subs(op(data), p_1), list)[1..3];
  else
    error("invalid point vector detected.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp := subs(op(data), p_2);
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := convert(subs(op(data), p_2), list)[1..3];
  else
    error("invalid point vector detected.");
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
  color::string               := TrussMe_FEM:-m_ElementColor
  }, $)::function;

  description "Plot the element from diplacements <d_1> and <d_2> given a list "
    "or set of data for substitution <data> and a display color <color>.";

  local d_1_tmp, d_2_tmp, p_1_tmp, p_2_tmp, p1p2, p1p2_unit, x, y, z,
    deformed_element;

  if type(d_1, list) and evalb(nops(d_1) = 6) then
    d_1_tmp := convert(d_1, Vector);
  elif type(d_1, Vector) and evalb(LinearAlgebra:-Dimension(d_1) = 6) then
    d_1_tmp := d_1;
  else
    error("invalid displacement vector detected.");
  end if;

  if type(d_2, list) and evalb(nops(d_2) = 6) then
    d_2_tmp := convert(d_2, Vector);
  elif type(d_2, Vector) and evalb(LinearAlgebra:-Dimension(d_2) = 6) then
    d_2_tmp := d_2;
  else
    error("invalid displacement vector detected.");
  end if;

  if type(p_1, list) and evalb(nops(p_1) = 3) then
    p_1_tmp := <convert(p_1, Vector), 1>;
  elif type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 4) then
    p_1_tmp := p_1;
  else
    error("invalid point vector detected.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp := <convert(p_2, Vector), 1>;
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := p_2;
  else
    error("invalid point vector detected.");
  end if;

  p1p2 := p_2_tmp - p_1_tmp;
  p1p2_unit := try abs(p1p2 /~ TrussMe_FEM:-Norm2(p1p2[1..3])) catch: [1,0,0] end try;

  # Calculate intermediate points through interpolation
  x := (d_1_tmp[1]*(1-xi) + d_2_tmp[1]*xi) * p1p2_unit[1] +
       (d_1_tmp[1]*(1-3*xi^2+2*xi^3) + d_2_tmp[1]*(3*xi^2-2*xi^3)) *
       TrussMe_FEM:-Norm2(p1p2_unit[2..3]) +
       (d_1_tmp[6]*(xi-xi^2)*(1-xi) + d_2_tmp[6]*(xi^2-xi^3)*xi) *
       TrussMe_FEM:-Norm2(p1p2_unit[2..3]) -
       (d_1_tmp[5]*(xi-xi^2)*(1-xi) + d_2_tmp[5]*(xi^2-xi^3)*xi) *
       TrussMe_FEM:-Norm2(p1p2_unit[2..3]);
  y := (d_1_tmp[2]*(1-xi) + d_2_tmp[2]*xi) * p1p2_unit[2] +
       (d_1_tmp[2]*(1-3*xi^2+2*xi^3) + d_2_tmp[2]*(3*xi^2-2*xi^3)) *
       TrussMe_FEM:-Norm2(p1p2_unit[[1,3]]) -
       (d_1_tmp[6]*(xi-xi^2)*(1-xi) + d_2_tmp[6]*(xi^2-xi^3)*xi) *
       TrussMe_FEM:-Norm2(p1p2_unit[[1,3]]) +
       (d_1_tmp[4]*(xi-xi^2)*(1-xi) + d_2_tmp[4]*(xi^2-xi^3)*xi) *
       TrussMe_FEM:-Norm2(p1p2_unit[[1,3]]);
  z := (d_1_tmp[3]*(1-xi) + d_2_tmp[3]*xi) * p1p2_unit[3] +
       (d_1_tmp[3]*(1-3*xi^2+2*xi^3) + d_2_tmp[3]*(3*xi^2-2*xi^3)) *
       TrussMe_FEM:-Norm2(p1p2_unit[1..2]) -
       (d_1_tmp[4]*(xi-xi^2)*(1-xi) + d_2_tmp[4]*(xi^2-xi^3)*xi) *
       TrussMe_FEM:-Norm2(p1p2_unit[1..2]) +
       (d_1_tmp[5]*(xi-xi^2)*(1-xi) + d_2_tmp[5]*(xi^2-xi^3)*(xi)) *
       TrussMe_FEM:-Norm2(p1p2_unit[1..2]);

  deformed_element := scaling *~ <x, y, z, 0> + p_1_tmp + p1p2 *~ xi;

  return plots:-display(
    plots:-spacecurve(
      subs(op(data), convert(deformed_element, list)[1..3]), xi = 0..1,
      parse("thickness") = 6
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
  color::string               := TrussMe_FEM:-m_LoadColor,
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
    error("invalid point vector detected.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp :=  subs(op(data), p_2);
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := convert( subs(op(data), p_2), list)[1..3];
  else
    error("invalid point vector detected.");
  end if;

  wb := max(0.02, 0.05 * scaling);
  wh := max(0.04, 0.10 * scaling);
  hh := min(0.1, try
      (0.1 * scaling)/TrussMe_FEM:-Norm2(p_2_tmp - p_1_tmp) catch: infinity end try
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
    "<nodes_frame_scaling, elements_frame_scaling>, and a loads scaling "
    "factor <scaling>.";

  local nodes, elements, loads, i, j, k, p_1, p_2, nodes_frame_scaling,
    elements_frame_scaling, disp_nodes_frames, disp_elements_frames, disp_nodes,
    disp_elements, disp_loads;

  # Extract the nodes, elements and loads from <fem> object
  nodes    := fem["nodes"];
  elements := fem["elements"];
  loads    := fem["loads"];

  # Retrieve frame scaling factors
  if type(frame_scaling, numeric) then
    nodes_frame_scaling := frame_scaling;
    elements_frame_scaling := frame_scaling;
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 1) then
    nodes_frame_scaling := frame_scaling[1];
    elements_frame_scaling := frame_scaling[1];
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 2) then
    nodes_frame_scaling := frame_scaling[1];
    elements_frame_scaling := frame_scaling[2];
  else
    error("<frame_scaling> expects a list of 1 or 2 elements.");
  end if;

  # Draw nodes frames
  if evalb(nodes_frame_scaling > 0) then
    disp_nodes_frames := [seq(i, i = 1..nops(nodes))];
    for i from 1 to nops(nodes) do
      disp_nodes_frames[i] := TrussMe_FEM:-DrawFrame(
        nodes[i]["frame"], parse("scaling") = nodes_frame_scaling
      );
    end do;
  else
    disp_nodes_frames := [];
  end if;

  # Draw elements frames
  if evalb(elements_frame_scaling > 0) then
    disp_elements_frames := [seq(i, i = 1..nops(elements))];
    for i from 1 to nops(elements) do
      disp_elements_frames[i] := TrussMe_FEM:-DrawFrame(
        elements[i]["frame"], parse("scaling") = elements_frame_scaling
      );
    end do;
  else
    disp_elements_frames := [];
  end if;

  # Plot the nodes
  disp_nodes := [seq(i, i = 1..nops(nodes))];
  for i from 1 to nops(nodes) do
    p_1 := <nodes[i]["coordinates"], 1>;
    disp_nodes[i] := TrussMe_FEM:-PlotNode(
      convert(p_1, Vector), parse("data")  = data,
      parse("color") = `if`(
        type(nodes[i], SUPPORT), TrussMe_FEM:-m_SupportColor, TrussMe_FEM:-m_NodeColor
      ),
      parse("token") = `if`(
        type(nodes[i], SUPPORT), TrussMe_FEM:-m_SupportToken, TrussMe_FEM:-m_NodeToken
      )
    );
  end do;

  # Plot the elements
  disp_elements := [seq(i, i = 1..nops(elements))];
  for i from 1 to nops(elements) do
    j := TrussMe_FEM:-GetObjById(nodes, elements[i]["node_1"], parse("position") = true);
    k := TrussMe_FEM:-GetObjById(nodes, elements[i]["node_2"], parse("position") = true);
    p_1 := <nodes[j]["coordinates"], 1>;
    p_2 := <nodes[k]["coordinates"], 1>;
    disp_elements[i] := TrussMe_FEM:-PlotElement(
      convert(p_1, Vector), convert(p_2, Vector), parse("data") = data,
      parse("color") = TrussMe_FEM:-m_ElementColor
    );
  end do;

  # Plot the loads
  if evalb(load_scaling > 0) then
    disp_loads := [seq(i, i = 1..2*nops(loads))];
    for i from 1 to nops(loads) do
      j := TrussMe_FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
      p_2 := <nodes[j]["coordinates"], 1>;

      # Plot forces
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      end if;
      disp_loads[2*i-1] := TrussMe_FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe_FEM:-m_ForceColor
      );

      # Plot moments
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      end if;
      disp_loads[2*i] := TrussMe_FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe_FEM:-m_MomentColor
      );
    end do;
  else
    disp_loads := [];
  end if;

  # Plot the structure
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
    "<nodes_frame_scaling, elements_frame_scaling>, a loads scaling factor "
    "<load_scaling>, and a deformation magnification factor "
    "<deformation_scaling>.";

  local nodes, elements, loads, i, j, k, p_1, p_2, d_1, d_2, nodes_frame_scaling,
    elements_frame_scaling, disp_nodes_frames, disp_elements_frames, disp_nodes,
    disp_elements, disp_loads, data_tmp;

  # Extract the nodes, elements and loads from <fem> object
  nodes    := fem["nodes"];
  elements := fem["elements"];
  loads    := fem["loads"];

  # Add data evaluated and reversed veils to data substitution
  data_tmp := [op(subs(op(data), ListTools:-Reverse(fem["veils"]))), op(data)];

  # Retrieve frame scaling factors
  if type(frame_scaling, numeric) then
    nodes_frame_scaling := frame_scaling;
    elements_frame_scaling := frame_scaling;
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 1) then
    nodes_frame_scaling := frame_scaling[1];
    elements_frame_scaling := frame_scaling[1];
  elif type(frame_scaling, list) and evalb(nops(frame_scaling) = 2) then
    nodes_frame_scaling := frame_scaling[1];
    elements_frame_scaling := frame_scaling[2];
  else
    error("<frame_scaling> expects a list of 1 or 2 elements.");
  end if;

  # Draw nodes frames
  if evalb(nodes_frame_scaling > 0) then
    disp_nodes_frames := [seq(i, i = 1..nops(nodes))];
    for i from 1 to nops(nodes) do
      # TODO: translate and rotate according to node deformations
      disp_nodes_frames[i] := TrussMe_FEM:-DrawFrame(
        nodes[i]["frame"], parse("scaling") = nodes_frame_scaling
      );
    end do;
  else
    disp_nodes_frames := [];
  end if;

  # Draw elements frames
  if evalb(elements_frame_scaling > 0) then
    disp_elements_frames := [seq(i, i = 1..nops(elements))];
    for i from 1 to nops(elements) do
      # TODO: translate and rotate according to node deformations
      disp_elements_frames[i] := TrussMe_FEM:-DrawFrame(
        elements[i]["frame"], parse("scaling") = elements_frame_scaling
      );
    end do;
  else
    disp_elements_frames := [];
  end if;

  # Plot the deformed nodes
  disp_nodes := [seq(i, i = 1..nops(nodes))];
  for i from 1 to nops(nodes) do
    p_1 := convert(
      <nodes[i]["coordinates"] +
        deformation_scaling *~ nodes[i]["frame"][1..3, 1..3].nodes[i]["output_displacements"][1..3], 1>,
      Vector
    );
    disp_nodes[i] := TrussMe_FEM:-PlotNode(
      p_1, parse("data")  = data_tmp,
      parse("color") = `if`(
        type(nodes[i], SUPPORT), TrussMe_FEM:-m_SupportColor, TrussMe_FEM:-m_NodeColor
      ),
      parse("token") = `if`(
        type(nodes[i], SUPPORT), TrussMe_FEM:-m_SupportToken, TrussMe_FEM:-m_NodeToken
      )
    );
  end do;

  # Plot the deformed elements
  disp_elements := [seq(i, i = 1..nops(elements))];
  for i from 1 to nops(elements) do
    j := TrussMe_FEM:-GetObjById(nodes, elements[i]["node_1"], parse("position") = true);
    k := TrussMe_FEM:-GetObjById(nodes, elements[i]["node_2"], parse("position") = true);
    if not interpolate then
      p_1 := <nodes[j]["coordinates"] +
        nodes[j]["frame"][1..3, 1..3].(deformation_scaling *~ nodes[j]["output_displacements"][1..3]), 1>;
      p_2 := <nodes[k]["coordinates"] +
        nodes[k]["frame"][1..3, 1..3].(deformation_scaling *~ nodes[k]["output_displacements"][1..3]), 1>;
      disp_elements[i] := TrussMe_FEM:-PlotElement(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data") = data_tmp, parse("color") = TrussMe_FEM:-m_ElementColor
      );
    else
      p_1 := <nodes[j]["coordinates"], 1>;
      p_2 := <nodes[k]["coordinates"], 1>;
      d_1 := <nodes[j]["frame"][1..3, 1..3].nodes[j]["output_displacements"][1..3],
        nodes[j]["frame"][1..3, 1..3].nodes[j]["output_displacements"][4..6]>;
      d_2 := <nodes[k]["frame"][1..3, 1..3].nodes[k]["output_displacements"][1..3],
        nodes[k]["frame"][1..3, 1..3].nodes[k]["output_displacements"][4..6]>;
      disp_elements[i] := TrussMe_FEM:-PlotDeformedElement(
        convert(p_1, Vector), convert(p_2, Vector),
        convert(d_1, Vector), convert(d_2, Vector),
        parse("data")    = data_tmp,
        parse("scaling") = deformation_scaling,
        parse("color")   = TrussMe_FEM:-m_ElementColor
      );
    end if;
  end do;

  # Plot the loads
  if evalb(load_scaling > 0) then
    disp_loads := [seq(i, i = 1..2*nops(loads))];
    for i from 1 to nops(loads) do
      j := TrussMe_FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
      p_2 := <nodes[j]["coordinates"] +
        deformation_scaling *~ nodes[j]["output_displacements"][1..3], 1>;

      # Plot forces
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
      end if;
      disp_loads[2*i-1] := TrussMe_FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data_tmp,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe_FEM:-m_ForceColor
      );

      # Plot moments
      if type(loads[i]["frame"], string) then
        p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      else
        p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
      end if;
      disp_loads[2*i] := TrussMe_FEM:-PlotLoad(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data")    = data_tmp,
        parse("scaling") = load_scaling,
        parse("color")   = TrussMe_FEM:-m_MomentColor
      );
    end do;
  else
    disp_loads := [];
  end if;

  # Plot the deformed structure
  return plots:-display(
    [op(disp_nodes_frames), op(disp_elements_frames),
     op(disp_nodes), op(disp_elements), op(disp_loads)],
    parse("axes")    = boxed,
    parse("scaling") = constrained,
    parse("labels")  = ['x', 'y', 'z']
  );
end proc: # PlotDeformedStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
