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
    return TrussMe_FEM:-m_NodeColor;
  elif type(obj, NODE) and not type(obj, SUPPORT) then
    return TrussMe_FEM:-m_SupportColor;
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

export PlotLoad := proc(
  p_1::POINT,
  p_2::POINT,
  {
  data::{list(`=`), set(`=`)} := [],
  color::string               := TrussMe_FEM:-m_LoadColor
  }, $)::function;

  description "Plot the load arrow from point <p_1> and <p_2> given a list or "
    "set of data for substitution <data> and a display color <color>.";

  # FIXME: width of the body, width and height of the head to be tuned
  return plots:-display(
    plottools:-arrow(
      convert(subs(op(data), p_1[1..3]), list),
      convert(subs(op(data), p_2[1..3]), list),
      0.08, 0.16, 0.04, cylindrical_arrow
    ),
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
    j := TrussMe_FEM:-GetObjById(
      nodes, elements[i]["node_1"], parse("position") = true
    );
    k := TrussMe_FEM:-GetObjById(
      nodes, elements[i]["node_2"], parse("position") = true
    );
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
  magnify::numeric            := 1,
  scaling::numeric            := 1,
  data::{list(`=`), set(`=`)} := []
  }, $)::{function, list(function)};

  description "Plot the undeformed structure made by <nodes>, <elements> and "
    "<loads> given a list or set of substitution data <data>, a loads scaling "
    "factor <scaling>, and a deformation magnification factor <magnify>.";


end proc: # PlotDeformedStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
