from bokeh.models import widgets as wid
from bokeh.plotting import figure
from bokeh.models import GraphRenderer, StaticLayoutProvider, Circle, HoverTool, BoxZoomTool, ResetTool, MultiLine
from bokeh.palettes import Spectral8, Inferno8
from bokeh.models.graphs import NodesAndLinkedEdges, EdgesAndLinkedNodes, from_networkx
from numpy import pi as pi
import networkx as nx
from numpy import ndarray, array
import bokeh.layouts as layout
from math import cos, sin
from functools import partial
from numpy import asarray


### Draw quadratic bezier paths
def bezier(start, end, control, steps):
    return [(1-s)**2*start + 2*(1-s)*s*control + s**2*end for s in steps]

def fullPlot(edges, names, cutoff):
    """
    fullPlot creates the fully connected and interactive graph plot given the parameters
    
    :param edges: contains all edges which are above the cutoff
    :param names: the names of all of the nodes
    :param cutoff: the cutoff which was used to assemble this particular plot
    :return p: the finished and interactive plot
    """
    
    num = len(names)
    indices = list(range(num))
   
    #Nodes are colored based on how many outgoing connectiosn they have normalized by the max number of outgoing nodes:            
    p = figure(title = "Full Learned Graph", x_range = (-1.1,1.1), y_range = (-1.1,1.1))
    p.plot_width = 700
    p.plot_height = 700
    p.title.text_font_size = "15px"
    
    G = nx.Graph()
    G.add_nodes_from(indices)
    edge_attrs = {}
    undirected, directed = "black", "red"
    outgoing = sum(edges)
    outgoing = ndarray.flatten(array(outgoing))
    
    for j in range(num):
        G.node[j]['name'] = names[j]
        G.node[j]['numOut'] = str(outgoing[j])
        for k in range(num):
            if edges[j,k] > cutoff:
                G.add_edge(j,k)
                edge_attrs[(j,k)] = directed
                if edges[k,j] > cutoff:
                    edge_attrs[(j,k)] = undirected
            
    #nx.set_edge_attributes(G, edge_attrs, 'edgeColor')
    
    norm = max(outgoing) #the normalizing constant to help with colors
    for i in range(num):
        if outgoing[i] == 0:
            G.node[i]['colors'] = Spectral8[0]
            G.node[i]['sColors'] = Inferno8[0]
        else:
            G.node[i]['colors'] = Spectral8[int(len(Spectral8)*((outgoing[i] - .001)/norm))]
            G.node[i]['sColors'] = Inferno8[int(len(Inferno8)*((outgoing[i] - .001)/norm))]
            
    node_hover_tool = HoverTool(tooltips=[("Name", "@name"), ("Num Outgoing", "@numOut")])
    p.add_tools(node_hover_tool, BoxZoomTool(), ResetTool())
    
    graph_renderer = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))

    graph_renderer.node_renderer.glyph = Circle(size = 14, fill_color = 'colors', name = 'name', 
                                                tags = ['numOut', nx.get_node_attributes(G,'numOut')])
    graph_renderer.node_renderer.selection_glyph = Circle(size = 16, fill_color = 'sColors')
    graph_renderer.node_renderer.hover_glyph = Circle(size = 16, fill_color = 'sColors')

    graph_renderer.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=4)
    graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color = '#CCCCCC', line_width=6)
    graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color = '#CCCCCC', line_width = 6)
    #graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color = 'edgeColor', line_width=6)
    #graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color = 'edgeColor', line_width = 6)
    
    graph_renderer.selection_policy = NodesAndLinkedEdges()
    graph_renderer.inspection_policy = EdgesAndLinkedNodes()
    
    p.renderers.append(graph_renderer)
    
    return(p)

def fullUpdate(pDict):
    """
    When the user updates the plot (by clicking on the button) this method will refresh the full plot)
    
    :param pDict: A dictionary containing all of the relevent information needed for updating - it is 
        automatically assembled by the code which calls it.  
    """
    
    tab = pDict['tab']; edges = pDict['edges']; names = pDict['names']; sel = pDict['selBox'].value
    cutoff = pDict['cutoff'].value
    
    inds = []
    for i in range(len(names)):
        if names[i] in set(sel):
            inds.append(i)
    xs = asarray([[i] for i in inds])
    ys = asarray([i for i in inds])
    newEdges = edges[xs, ys]
    newNames = [names[i] for i in inds]
    
    fPlot = fullPlot(newEdges, newNames, cutoff)
    
    org = tab.child
    org.children[1] = fPlot #update the plot in the tab  
    tab.child = org 

def fullHelp(hDict):
    """
    fullHelp updates a text box if the user selects one of the FAQ and clicks the help button
    
    :param hDict: contains all the parameters (such as which question was selected) that are 
        needed for updating the FAQ box.
    """
    tab = hDict['tab']; hReq = hDict['selected'].value; opts = hDict['options']
    
    if(hReq == opts[0]):
        tab.child.children[5].children[4] = wid.Div(text = """<p style = 
        "text-align:center;font-family:verdana"> 
        Please selecte a FAQ from the drop down to recieve help on!<p>""")
        
    if(hReq == opts[1]):
         tab.child.children[5].children[4] = wid.Div(text = """<p style = 
         "text-align:center;font-family:verdana"> 
         This plot shows the evolution of the probability that the 
         feature (edge) you selected is in the learned graph. As one 
         samples from more graphs within the posterior, the estimate 
         (as well as its confidence interval) evolve. <p>""")
    

def fullTab(edges, names):
    """
    fullTab assembles the tab structure for the fullGraph tab
    :param edges: a matrix with all edge likelihoods learned from minimal IMAP MCMC
    :param names: the names of all nodes learned
    :return tab2: the tab object containing all objects in the edgeMat tab
    """
    cutoff = wid.Slider(start = 0, end = 1, value = .5, step=.01, title="Cutoff")
    fullGraph = fullPlot(edges, names, cutoff.value)
    listBox = wid.MultiSelect(title = "Nodes to Show:", options = list(names), size = 20,
                              value = list(names))
    replotButton = wid.Button(label = "Replot", button_type = "success")
    instructions = wid.Div(text = """<p style = "text-align:center;font-family:verdana"> 
        Here is a representation of the full learned graph.  \n 
        To reduce the number of nodes shown, select your \n
        desired nodes from the multiselect tool.<p>""")
    sp0 = wid.Div(width = 60); sp3 = wid.Div(width = 100)
    sp1 = wid.Div(height = 40); sp2 = wid.Div(height = 10) #spacers to make the tab more visually appealing
    sp4 = wid.Div(width = 40); sp5 = wid.Div(height = 50)
    sp6 = wid.Div(height = 50)
    
    helpButton = wid.Button(label = "Help?")
    helpOptions = ["", "How do I read this plot?"]
    helpDrop = wid.Select(title = "FAQ:", value = helpOptions[0], options = helpOptions)
    helpRead = wid.Div(width = 40)
    
    tab3 = wid.Panel(child = layout.row(sp0, fullGraph, sp3, layout.column(sp1, instructions, sp2, cutoff, listBox,
                        replotButton),sp4, layout.column(sp5, helpButton, helpDrop, sp6, helpRead)), title = "Full Graph")
    
    passDict = {'tab': tab3, 'edges': edges, 'names': names, 'selBox': listBox, 'cutoff': cutoff}
    replotButton.on_click(partial(fullUpdate, pDict = passDict))
    helpDict = {'tab': tab3, 'selected': helpDrop, 'options': helpOptions}
    helpButton.on_click(partial(fullHelp, hDict = helpDict))

    return(tab3)