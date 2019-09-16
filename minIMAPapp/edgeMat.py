from bokeh import plotting as plt
from bokeh.models import ColorBar, LinearColorMapper
from bokeh.palettes import Spectral8 
from bokeh.models import widgets as wid
import bokeh.layouts as layout
from numpy import pi, asarray
from functools import partial

def edgePlot(gData, names, cutoff):
    """
    edgePlot generates the heat map which corresponds to all the edges present in the DAG.
    :param gData: takes in the data dict which contains all the values used for the edge informtaion and names etc
    :param names: the names of the various nodes being looked at 
    
    :return p: the bokeh figure associated with the edges
    """
    
    p = plt.figure(title = "Features Likely Present in DAG with Frequency above " + str(cutoff),
                   x_axis_location = "above", tools = "hover, save, box_zoom, wheel_zoom, pan, reset",
                   x_range = list(reversed(names)), y_range = names,
                   tooltips = [('Dependence', '@yname @dep @xname'), ('Prob in graph:', '@pval')])
    
    #set a few things about the plot appearance:
    p.plot_width = 700
    p.plot_height = 700
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi/3
    p.title.text_font_size = "15px"
    p.title.align = 'center'
    
    #create the little rectangles which create the heatmap:
    p.rect('xname', 'yname', 0.9, 0.9, source = gData,
           color = 'colors', alpha = 'alphas', line_color=None,
           hover_line_color='black', hover_color = 'colors')
    cMap = LinearColorMapper(palette = Spectral8, low = cutoff, high = 1)
    colorbar = ColorBar(color_mapper = cMap, location = (0,0), label_standoff = 10)
    p.add_layout(colorbar, 'right')
    
    return(p)


def edgePackage(edges, names, cutoff):
    """
    edgePackage packages the edge information into a form which can be read and utilized by edgeMat viewer.       
    
    :param edges: a matrix containing all edge likelihoods learned from minIMAP MCMC
    :param names: the names of the nodes
    :return graphDict: a dictionary with components for the graphical program.  Contains info about the p values of various 
        links, as well as info regarding their color, dependency etc etc.
    """
    
    colors = []; xNames = []; yNames = []; alphas = []; dep = []; pval = []
    size = len(names) #the number of nodes in the graph
    
    for i in range(size):
        for j in range(size):
        
            xNames.append(names[i])
            yNames.append(names[j])           
            
            if edges[i,j] > cutoff:
                if edges[j,i] > cutoff: #ie if the directionality of the edge cannot be determined:
                    dep.append('share a connection ')
                    colors.append(Spectral8[int(len(Spectral8)*((1-edges[i,j])/(1-cutoff)))])
                    alphas.append(.5)
                    pval.append(str(edges[i,j]))
                else: # likely a directed edge                 
                    dep.append('likely depends on ')
                    colors.append(Spectral8[int(len(Spectral8)*((1-edges[i,j])/(1-cutoff)))])
                    alphas.append(1)
                    pval.append(str(edges[i,j]))
            else: 
                colors.append('lightgray')
                dep.append('likely doesn\'t depend on')
                alphas.append(.5)
                pval.append(' <= ' + str(cutoff))

    graphData = dict(xname = xNames, yname = yNames, colors = colors, alphas = alphas, pval = pval, dep = dep)
    return(graphData)

def edgeUpdate(pDict):
    """
    When the user updates the plot (by clicking on the button) this method will refresh the edge plot)
    
    :param pDict: A dictionary containing all of the relevent information needed for updating - it is 
        automatically assembled by the code which calls it.  
    """
    
    tab = pDict['tab']; edges = pDict['edges'].copy(); names = pDict['names'].copy(); valSel = pDict['selected'].value
    cutoff = pDict['cutoff'].value
    inds = []
    for i in range(len(names)):
        if names[i] in set(valSel):
            inds.append(i)
    xs = asarray([[i] for i in inds])
    ys = asarray([i for i in inds])
    edges = edges[xs, ys]
    names = [names[i] for i in inds]
    newPack = edgePackage(edges, names, cutoff)
    ePlot = edgePlot(newPack, names, cutoff)
    org = tab.child
    org.children[1] = ePlot #update the plot in the tab  
    tab.child = org

def edgeHelp(hDict):
    """
    edgeHelp updates a text box if the user selects one of the FAQ and clicks the help button
    
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
         This plot represents an edge connection matrix.  
         To read it, choose a square of interest.  
         If colored in, that means the probability that there is an edge 
         between the variable representing the column that the square is in
         (to read its name go to the top) and the variable of its row
         (whose name is to the left of the matrix) is above the cutoff you have
         set.  If the square is instead light gray, then the probability of this
         feature does not surpass the cutoff.  In practice it might be easier to
         scroll over the square to have these values read off.<p>""")
    

def edgeTab(edges, names):
    """
    edgeTab assembles the tab structure for the EdgeMat tab
    
    :param edges: a matrix with all edge likelihoods learned from minimal IMAP MCMC
    :param names: the names of all nodes learned
    :return tab1: the tab object containing all objects in the edgeMat tab
    """
    cutoff = wid.Slider(start = 0, end = 1, value = .5, step=.01, title="Cutoff")
    gData = edgePackage(edges, names, cutoff.value)
    eGraph = edgePlot(gData, names, cutoff.value)
    listBox = wid.MultiSelect(title = "Nodes to Show:", options = list(names), size = 20,
                              value = list(names))
    sp0 = wid.Div(width = 60); sp3 = wid.Div(width = 100)
    sp1 = wid.Div(height = 40); sp2 = wid.Div(height = 10); #spacers to make the tab more visually appealing
    sp4 = wid.Div(width = 100); sp5 = wid.Div(height = 40)
    sp6 = wid.Div(height = 60)
    replotButton = wid.Button(label = "Replot", button_type = "success")
    helpButton = wid.Button(label = "Help?")
    helpOptions = ["", "How do I read this plot?"]
    helpDrop = wid.Select(title = "FAQ:", value = helpOptions[0], options = helpOptions)
    helpRead = wid.Div(width = 40)
    instructions = wid.Div(text = """<p style = "text-align:center;font-family:verdana"> 
        Colored boxes represent edges from the column variable to the row variable.  
        To view a smaller subsample of the nodes, select them using the tool below.<p>""")
    
    tab1 = wid.Panel(child = layout.row(sp0, eGraph, sp3, layout.column(sp1, instructions, sp2, cutoff, listBox, replotButton), 
                                       sp4, layout.column(sp5, helpButton, helpDrop, sp6, helpRead)),
                 title = "Connections Matrix")
    passDict = {'tab': tab1, 'selected': listBox, 'names': names, 'edges': edges, 'cutoff': cutoff}
    helpDict = {'tab': tab1, 'selected': helpDrop, 'options': helpOptions}
    replotButton.on_click(partial(edgeUpdate, pDict = passDict))
    helpButton.on_click(partial(edgeHelp, hDict = helpDict))

    return(tab1)