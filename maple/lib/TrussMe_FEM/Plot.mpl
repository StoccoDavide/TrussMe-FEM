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
    p_tmp := p;
  elif type(p, Vector) and evalb(LinearAlgebra:-Dimension(p) = 4) then
    p_tmp := convert(p, list)[1..3];
  else
    error("invalid point vector detected.");
  end if;

  return plots:-display(
    plottools:-point(
      subs(op(data), p_tmp),
      parse("symbol")     = token,
      parse("symbolsize") = 15
    ),
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
    p_1_tmp := p_1;
  elif type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 4) then
    p_1_tmp := convert(p_1, list)[1..3];
  else
    error("invalid point vector detected.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp := p_2;
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := convert(p_2, list)[1..3];
  else
    error("invalid point vector detected.");
  end if;

  return plots:-display(
    plottools:-line(
      subs(op(data), p_1_tmp),
      subs(op(data), p_2_tmp),
      parse("thickness") = 6
    ),
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
  magnify::nonnegative        := 1.0,
  data::{list(`=`), set(`=`)} := [],
  color::string               := TrussMe_FEM:-m_ElementColor
  }, $)::function;

  description "Plot the element from diplacements <d_1> and <d_2> given a list "
    "or set of data for substitution <data> and a display color <color>.";

  local d_1_tmp, d_2_tmp, p_1_tmp, p_2_tmp, p1p2, p1p2_unit, x, y, z, deformed_element;

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
  p1p2_unit := abs(p1p2 /~ TrussMe_FEM:-Norm2(p1p2[1..3]));

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

  deformed_element := magnify *~ <x, y, z, 0> + p_1_tmp + p1p2 *~ xi;

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
    p_1_tmp := p_1;
  elif type(p_1, Vector) and evalb(LinearAlgebra:-Dimension(p_1) = 4) then
    p_1_tmp := convert(p_1, list)[1..3];
  else
    error("invalid point vector detected.");
  end if;

  if type(p_2, list) and evalb(nops(p_2) = 3) then
    p_2_tmp := p_2;
  elif type(p_2, Vector) and evalb(LinearAlgebra:-Dimension(p_2) = 4) then
    p_2_tmp := convert(p_2, list)[1..3];
  else
    error("invalid point vector detected.");
  end if;

  p_1_tmp := subs(op(data), p_1_tmp);
  p_2_tmp := subs(op(data), p_2_tmp);
  wb := max(0.02, 0.05 * scaling);
  wh := max(0.04, 0.10 * scaling);
  hh := min(0.1, try (0.1 * scaling)/TrussMe_FEM:-Norm2(p_2_tmp - p_1_tmp) catch: infinity end try);
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
  load_scaling::numeric            := 1,
  data::{list(`=`), set(`=`)} := []
  }, $)::{function, list(function)};

  description "Plot the undeformed structure made by <nodes>, <elements> and "
    "<loads> given a list or set of substitution data <data>, and a loads "
    "scaling factor <scaling>.";

  local i, j, k, p_1, p_2, disp_nodes, disp_elements, disp_loads;

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
      parse("data") = data, parse("color") = TrussMe_FEM:-m_ForceColor,
      parse("scaling") = load_scaling
    );

    # Plot moments
    if type(loads[i]["frame"], string) then
      p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
    else
      p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
    end if;
    disp_loads[2*i] := TrussMe_FEM:-PlotLoad(
      convert(p_1, Vector), convert(p_2, Vector),
      parse("data") = data, parse("color") = TrussMe_FEM:-m_MomentColor,
      parse("scaling") = load_scaling
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
  interpolate::boolean        := true,
  magnify::numeric            := 1,
  load_scaling::numeric       := 1,
  data::{list(`=`), set(`=`)} := []
  }, $)::{function, list(function)};

  description "Plot the undeformed structure made by <nodes>, <elements> and "
    "<loads> given a list or set of substitution data <data>, a loads scaling "
    "factor <scaling>, and a deformation magnification factor <magnify>.";

  local i, j, k, p_1, p_2, d_1, d_2, disp_nodes, disp_elements, disp_loads;

  # Plot the deformed nodes
  disp_nodes := [seq(i, i = 1..nops(nodes))];
  for i from 1 to nops(nodes) do
    p_1 := convert(
      <convert(nodes[i]["coordinates"], Matrix) +
        magnify *~ nodes[i]["frame"][1..3, 1..3].nodes[i]["output_displacements"][1..3], 1>,
      Vector
    );
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
    if not interpolate then
      p_1 := <convert(nodes[j]["coordinates"], Matrix) +
        magnify *~ nodes[j]["frame"][1..3, 1..3].nodes[j]["output_displacements"][1..3], 1>;
      p_2 := <convert(nodes[k]["coordinates"], Matrix) +
        magnify *~ nodes[k]["frame"][1..3, 1..3].nodes[k]["output_displacements"][1..3], 1>;
      disp_elements[i] := TrussMe_FEM:-PlotElement(
        convert(p_1, Vector), convert(p_2, Vector),
        parse("data") = data, parse("color") = TrussMe_FEM:-m_ElementColor
      );
    else
      p_1 := <nodes[j]["coordinates"], 1>;
      p_2 := <nodes[k]["coordinates"], 1>;
      d_1 := <nodes[j]["frame"][1..3, 1..3].nodes[j]["output_displacements"][1..3],
        nodes[j]["frame"][1..3, 1..3].nodes[j]["output_displacements"][4..6]>;
      d_2 := <nodes[k]["frame"][1..3, 1..3].nodes[k]["output_displacements"][1..3],
        nodes[k]["frame"][1..3, 1..3].nodes[k]["output_displacements"][4..6]>;
      disp_elements[i] := TrussMe_FEM:-PlotDeformedElement(
        convert(p_1, Vector), convert(p_2, Vector), convert(d_1, Vector), convert(d_2, Vector),
        parse("data") = data, parse("color") = TrussMe_FEM:-m_ElementColor,
        parse("magnify") = magnify
      );
    end if;
  end do;

  # Plot the loads
  disp_loads := [seq(i, i = 1..2*nops(loads))];
  for i from 1 to nops(loads) do
    j := TrussMe_FEM:-GetObjById(nodes, loads[i]["node"], parse("position") = true);
    p_2 := <convert(nodes[j]["coordinates"], Matrix) +
      magnify *~ nodes[j]["output_displacements"][1..3], 1>;

    # Plot forces
    if type(loads[i]["frame"], string) then
      p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
    else
      p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][1..3], 0>;
    end if;
    disp_loads[2*i-1] := TrussMe_FEM:-PlotLoad(
      convert(p_1, Vector), convert(p_2, Vector),
      parse("data") = data, parse("color") = TrussMe_FEM:-m_ForceColor,
      parse("scaling") = load_scaling
    );

    # Plot moments
    if type(loads[i]["frame"], string) then
      p_1 := p_2 - nodes[j]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
    else
      p_1 := p_2 - loads[i]["frame"].<load_scaling * loads[i]["components"][4..6], 0>;
    end if;
    disp_loads[2*i] := TrussMe_FEM:-PlotLoad(
      convert(p_1, Vector), convert(p_2, Vector),
      parse("data") = data, parse("color") = TrussMe_FEM:-m_MomentColor,
      parse("scaling") = load_scaling
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
