from numpy import mod, floor, asarray
from bokeh import plotting as plt
from numpy import ndarray, array
import bokeh.layouts as layout
from bokeh.models import widgets as wid
from functools import partial

def learnUpdate(pDict):
    
    tab = pDict['tab']; evolve = pDict['evolve']; startV = pDict['start'].value; endV = pDict['end'].value
    names = pDict['names']
    start = 0; end = 0
    for i in range(len(names)):
        if names[i] == startV:
            start = i
        if names[i] == endV:
            end = i
    ind = start*len(names) + end
    edge = ndarray.flatten(array(evolve[ind]))
    tPlot = learnPlot(startV, endV, edge)
    
    credLab = wid.Div(text = """<p style = 
         "text-align:center;font-family:verdana"> 
         Probability Edge is On: """ + str(edge[-1]) + """ <p>""")
    
    org = tab.child
    org.children[3] = tPlot #update the plot in the tab 
    org.children[1].children[1] = credLab
    tab.child = org

def learnPlot(sName, eName, edge):
    
    p = plt.figure(title = "Evolution of prob of edge from " + sName + " to  " + eName)
    p.plot_width = 700
    p.plot_height = 700
    p.title.text_font_size = "15px"
    confY, confX = confInt(ndarray.flatten(array(edge)), 5)
    p.patch(x = confX, y = confY, color = '#99d8c9')
    p.line(asarray([i for i in range(len(edge))]), ndarray.flatten(array(edge)))
    p.yaxis.axis_label = "Probability of Edge from "  + sName + " to  " + eName
    p.xaxis.axis_label = "Number of Samples from the Posterior"
    
    return(p)

def learnHelp(hDict):
    
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
    

def confInt(pEst, stepSz):
    """
    confInt will return a series of x and y points to define a polygon patch represting the confidence interval over
    the estimate for any given edge.
    """
    k = 0
    yBot = []; xBot = [] 
    yTop = []; xTop = [] #lists to put polygon intergers in
    while (k*stepSz + 10) < len(pEst):
        i = k*stepSz + 10
        sd = (pEst[i]*(1 - pEst[i])/(i))**(1/2)
        yTop.append(pEst[i] + sd)
        xTop.append(i)
        yBot.insert(0, pEst[i] - sd)
        xBot.insert(0, i)
        k += 1
    
    return(yTop + yBot, xTop + xBot)
    
    
def learnTab(evolve, names):
    """
    learnTab is a diagnostics tab which displays the traceplots of the edges.  It shows the evolution in estimation of the
        probability of each feature over the course of learning.  Also shown are evolving CI intervals.
    
    :param evolve: the large matrix containing arrays of edge confidence values indexe by time
    :param names: the names of the nodes under consideration
    :return tab0: the tab containing all the information looped together    
    """
    start = 0; end = 1; ind = start*len(names) + end
    startSel = wid.Select(title = "Orgin of Edge:", value = names[start], options = names)
    endSel = wid.Select(title = "End of Edge:", value = names[end], options = names)
    rPlot = wid.Button(label = "Replot", button_type = "success")
    edge = ndarray.flatten(array(evolve[ind]))
    credLab = wid.Div(text = """<p style = 
         "text-align:center;font-family:verdana"> 
         Posterior Edge is On: """ + str(edge[-1]) + """ <p>""")
    
    sp0 = wid.Div(width = 60); sp1 = wid.Div(width = 100)
    sp2 = wid.Div(width = 60); sp3 = wid.Div(height = 130)
    sp4 = wid.Div(height = 100); sp5 = wid.Div(height = 60)
    sp6 = wid.Div(height = 60)
    
    helpButton = wid.Button(label = "Help?")
    helpOptions = ["", "How do I read this plot?"]
    helpDrop = wid.Select(title = "FAQ:", value = helpOptions[0], options = helpOptions)
    helpRead = wid.Div(width = 40)
    
    tracePlot = learnPlot(names[start], names[end], edge)
    
    spread = layout.row(sp0, layout.column(sp4,credLab, sp6, startSel, endSel, rPlot), sp1, tracePlot, sp2, 
                       layout.column(sp3, helpButton, helpDrop, sp5, helpRead))
    tab2 = wid.Panel(child = spread, title = "Evaluate Convergence")
    
    passDict = {'tab': tab2, 'evolve': evolve, 'start': startSel, 'end': endSel, 'names': names}
    rPlot.on_click(partial(learnUpdate, pDict = passDict))
    helpDict = {'tab': tab2, 'selected': helpDrop, 'options': helpOptions}
    helpButton.on_click(partial(learnHelp, hDict = helpDict))
    
    return(tab2)