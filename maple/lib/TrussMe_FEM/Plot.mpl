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
  nodes::NODES,
  elements::ELEMENTS,
  loads::LOADS := [],
  {
  id::boolean   := false,
  disp::boolean := true
  }, $)::function;

  description "Return the graph of the structure given a list of nodes "
    "<nodes> and elements <elements>.";

  local i, vertex_name, vertex_id, vertex_color, G;

  if m_VerboseMode then
    printf("TrussMe:-StructureGraph(...): checking structure connections... ");
  end if;

  # Retrieve the objects names, ids and colors
  vertex_name := [
    seq(nodes[i]["name"],    i = 1..nops(nodes)),
    seq(elements[i]["name"], i = 1..nops(elements)),
    seq(loads[i]["name"],    i = 1..nops(loads))
  ];
  vertex_id := [
    seq(nodes[i]["id"],    i = 1..nops(nodes)),
    seq(elements[i]["id"], i = 1..nops(elements)),
    seq(loads[i]["id"],    i = 1..nops(loads))
  ];
  vertex_color := [
    seq(TrussMe_FEM:-ObjectColor(nodes[i]), i = 1..nops(nodes)),
    seq(TrussMe_FEM:-ObjectColor(elements[i]), i = 1..nops(elements)),
    seq(TrussMe_FEM:-ObjectColor(loads[i]), i = 1..nops(loads))
  ];

  # Build the graph
  G := GraphTheory:-Graph(vertex_id);
  GraphTheory:-HighlightVertex(G, vertex_id, vertex_color);

  # Connect the nodes with the elements
  for i from 1 to nops(elements) do
    GraphTheory:-AddEdge(G, {elements[i]["id"], elements[i]["node_1"]});
    GraphTheory:-AddEdge(G, {elements[i]["id"], elements[i]["node_2"]});
  end do;

  # Connect the nodes with the loads
  for i from 1 to nops(loads) do
    GraphTheory:-AddEdge(G, {loads[i]["id"], loads[i]["node"]});
  end do;

  # Relabel with vertex names
  if not id then
    G := GraphTheory:-RelabelVertices(G, vertex_name);
  end if;

  if m_VerboseMode then
    printf("DONE\n");
  end if;

  # Plot the graph
  if disp then
    print(plots:-display(GraphTheory:-DrawGraph(G, layout = tree)));
  end if;

  # Return the graph
  return G;
end proc: # StructureGraph

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotNode := proc(
  p::POINT,
  {
  data::{list(`=`), set(`=`)} := [],
  token::symbol               := TrussMe_FEM:-m_NodeToken,
  color::string               := TrussMe_FEM:-m_NodeColor
  }, $)::function;

  description "Plot the node (or support) at point <p> given a list or set of "
    "data for substitution <data>, a display token <token> and a display color "
    "<color>.";

  return plots:-display(
    plottools:-point(
      convert(subs(op(data), p[1..3]), list),
      parse("symbol")     = token,
      parse("symbolsize") = 20
    ),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotNode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotElement := proc(
  p_1::POINT,
  p_2::POINT,
  {
  data::{list(`=`), set(`=`)} := [],
  color::string               := TrussMe_FEM:-m_ElementColor
  }, $)::function;

  description "Plot the element from point <p_1> and <p_2> given a list or set "
    "of data for substitution <data> and a display color <color>.";

  return plots:-display(
    plottools:-line(
      convert(subs(op(data), p_1[1..3]), list),
      convert(subs(op(data), p_2[1..3]), list),
      parse("thickness") = 6
    ),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedElement := proc(
  p_1::POINT,
  p_2::POINT,
  d_1::{list(algebraic), Vector(algebraic)},
  d_2::{list(algebraic), Vector(algebraic)},
  {
  frame::FRAME                := Matrix(4, shape = identity),
  divisions::nonnegint        := 10,
  magnify::nonnegative        := 1.0,
  data::{list(`=`), set(`=`)} := [],
  color::string               := TrussMe_FEM:-m_ElementColor
  }, $)::function;

  description "Plot the element from diplacements <d_1> and <d_2> given a list "
    "or set of data for substitution <data>, a display color <color> and the "
    "numer of element divisions <divisions>.";

  local d_1_tmp, d_2_tmp, x, y, z;

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

  # Add nodes positions to the displacement vectors
  d_1_tmp[1..3] := p_1 + magnify *~ d_1_tmp;

  # Calculate intermediate points through interpolation
  x := d_1_tmp[1]*(1-3*xi^2+2*xi^3) + d_1_tmp[4]*(xi-xi^2)*(1-xi) +
       d_2_tmp[1]*(3*xi^2-2*xi^3)   + d_2_tmp[4]*(xi^2-xi^3)*xi;
  y := d_1_tmp[2]*(1-3*xi^2+2*xi^3) + d_1_tmp[5]*(xi-xi^2)*(1-xi) +
       d_2_tmp[2]*(3*xi^2-2*xi^3)   + d_2_tmp[5]*(xi^2-xi^3)*xi;
  z := d_1_tmp[3]*(1-3*xi^2+2*xi^3) + d_1_tmp[6]*(xi-xi^2)*(1-xi) +
       d_2_tmp[3]*(3*xi^2-2*xi^3)   + d_2_tmp[6]*(xi^2-xi^3)*xi;

  return plots:-display(
    plots:-spacecurve(
      convert(subs(op(data), (frame.<x, y, z, 1>)[1..3]), list), xi = 0..1,
      parse("thickness") = 6
    ),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotDeformedElement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotLoad := proc(
  p_1::POINT,
  p_2::POINT,
  {
  data::{list(`=`), set(`=`)} := [],
  color::string               := TrussMe_FEM:-m_LoadColor
  }, $)::function;

  description "Plot the load arrow from point <p_1> and <p_2> given a list or "
    "set of data for substitution <data> and a display color <color>.";

  local p_1_tmp, p_2_tmp, wb, wh, hh;

  p_1_tmp := convert(subs(op(data), p_1[1..3]), list);
  p_2_tmp := convert(subs(op(data), p_2[1..3]), list);
  wb := 0.05;
  wh := 0.10;
  hh := max(0.1, 0.1/TrussMe_FEM:-Norm2(p_2_tmp - p_1_tmp));
  return plots:-display(
    plottools:-arrow(p_1_tmp, p_2_tmp, wb, wh, hh, cylindrical_arrow),
    parse("linestyle") = solid,
    parse("color")     = color,
    parse("scaling")   = constrained
  );
end proc: # PlotLoad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotStructure := proc(
  nodes::NODES       := [],
  elements::ELEMENTS := [],
  loads::LOADS       := [],
  {
  scaling::numeric            := 1,
  data::{list(`=`), set(`=`)} := []
  }, $)::{function, list(function)};

  description "Plot the undeformed structure made by <nodes>, <elements> and "
    "<loads> given a list or set of substitution data <data>, and a loads "
    "scaling factor <scaling>.";

  local i, j, k, p_1, p_2, disp_nodes, disp_elements, disp_loads;

  # Plot the nodes
  disp_nodes := [seq(i, i = 1..nops(nodes))];
  for i from 1 to nops(nodes) do
    p_1 := convert(nodes[i]["frame"].<nodes[i]["coordinates"], 1>, Vector);
    disp_nodes[i] := TrussMe_FEM:-PlotNode(
      p_1, parse("data")  = data,
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
    p_1 := convert(nodes[j]["frame"].<nodes[j]["coordinates"], 1>, Vector);
    p_2 := convert(nodes[k]["frame"].<nodes[k]["coordinates"], 1>, Vector);
    disp_elements[i] := TrussMe_FEM:-PlotElement(
      p_1, p_2, parse("data") = data, parse("color") = TrussMe_FEM:-m_ElementColor
    );
  end do;

  # Plot the loads
  disp_loads := [seq(i, i = 1..2*nops(loads))];
  for i from 1 to nops(loads) do
    j := TrussMe_FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
    p_2 := nodes[j]["frame"].<nodes[j]["coordinates"], 1>;

    # Plot forces
    if type(loads[i]["frame"], string) then
      p_1 := p_2 - nodes[j]["frame"].<scaling * loads[i]["components"][1..3], 0>;
    else
      p_1 := p_2 - loads[i]["frame"].<scaling * loads[i]["components"][1..3], 0>;
    end if;
    disp_loads[2*i-1] := TrussMe_FEM:-PlotLoad(
      convert(p_1, Vector), convert(p_2, Vector),
      parse("data") = data, parse("color") = TrussMe_FEM:-m_ForceColor
    );

    # Plot moments
    if type(loads[i]["frame"], string) then
      p_1 := p_2 - nodes[j]["frame"].<scaling * loads[i]["components"][4..6], 0>;
    else
      p_1 := p_2 - loads[i]["frame"].<scaling * loads[i]["components"][4..6], 0>;
    end if;
    disp_loads[2*i] := TrussMe_FEM:-PlotLoad(
      convert(p_1, Vector), convert(p_2, Vector),
      parse("data") = data, parse("color") = TrussMe_FEM:-m_MomentColor
    );
  end do;

  # Plot the structure
  return plots:-display(
    [op(disp_nodes), op(disp_elements), op(disp_loads)],
    parse("axes")    = boxed,
    parse("scaling") = constrained,
    parse("labels")  = ['x', 'y', 'z']
  );
end proc: # PlotStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedStructure := proc(
  nodes::NODES       := [],
  elements::ELEMENTS := [],
  loads::LOADS       := [],
  {
  divisions::nonnegint        := 0,
  magnify::numeric            := 1,
  scaling::numeric            := 1,
  data::{list(`=`), set(`=`)} := []
  }, $)::{function, list(function)};

  description "Plot the undeformed structure made by <nodes>, <elements> and "
    "<loads> given a list or set of substitution data <data>, a loads scaling "
    "factor <scaling>, and a deformation magnification factor <magnify>.";

  local i, j, k, p_1, p_2, disp_nodes, disp_elements, disp_loads;

  # Plot the deformed nodes
  disp_nodes := [seq(i, i = 1..nops(nodes))];
  for i from 1 to nops(nodes) do
    p_1 := convert(nodes[i]["frame"].<(
      convert(nodes[i]["coordinates"], Matrix) + magnify *~ nodes[i]["output_displacements"][1..3]
    ), 1>, Vector);
    disp_nodes[i] := TrussMe_FEM:-PlotNode(
      p_1, parse("data")  = data,
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
    p_1 := nodes[j]["frame"].<convert(nodes[j]["coordinates"], Matrix) +
      magnify *~ nodes[j]["output_displacements"][1..3], 1>;
    p_2 := nodes[k]["frame"].<convert(nodes[k]["coordinates"], Matrix) +
      magnify *~ nodes[k]["output_displacements"][1..3], 1>;
    if evalb(divisions = 0) then
      disp_elements[i] := TrussMe_FEM:-PlotElement(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data") = data, parse("color") = TrussMe_FEM:-m_ElementColor
      );
    else
      d_1 := magnify *~ nodes[j]["output_displacements"];
      d_2 := magnify *~ nodes[k]["output_displacements"];
      disp_elements[i] := TrussMe_FEM:-PlotDeformedElement(
        convert(p_1, Vector), convert(p_2, Vector),
        convert(d_1, Vector), convert(d_2, Vector),
        parse("data") = data, parse("color") = TrussMe_FEM:-m_ElementColor,
        parse("frame") = elements[i]["frame"], parse("divisions") = divisions,
        parse("magnify") = magnify
      );
    end if;
  end do;

  # Plot the loads
  disp_loads := [seq(i, i = 1..2*nops(loads))];
  for i from 1 to nops(loads) do
    j := TrussMe_FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
    p_2 := nodes[j]["frame"].<convert(nodes[j]["coordinates"], Matrix) +
      magnify *~ nodes[j]["output_displacements"][1..3], 1>;

    # Plot forces
    if type(loads[i]["frame"], string) then
      p_1 := p_2 - nodes[j]["frame"].<scaling * loads[i]["components"][1..3], 0>;
    else
      p_1 := p_2 - loads[i]["frame"].<scaling * loads[i]["components"][1..3], 0>;
    end if;
    disp_loads[2*i-1] := TrussMe_FEM:-PlotLoad(
      convert(p_1, Vector), convert(p_2, Vector),
      parse("data") = data, parse("color") = TrussMe_FEM:-m_ForceColor
    );

    # Plot moments
    if type(loads[i]["frame"], string) then
      p_1 := p_2 - nodes[j]["frame"].<scaling * loads[i]["components"][4..6], 0>;
    else
      p_1 := p_2 - loads[i]["frame"].<scaling * loads[i]["components"][4..6], 0>;
    end if;
    disp_loads[2*i] := TrussMe_FEM:-PlotLoad(
      convert(p_1, Vector), convert(p_2, Vector),
      parse("data") = data, parse("color") = TrussMe_FEM:-m_MomentColor
    );
  end do;

  # Plot the deformed structure
  return plots:-display(
    [op(disp_nodes), op(disp_elements), op(disp_loads)],
    parse("axes")    = boxed,
    parse("scaling") = constrained,
    parse("labels")  = ['x', 'y', 'z']
  );
end proc: # PlotDeformedStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
