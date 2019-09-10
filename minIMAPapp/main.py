from edgeMat import *
from learnPlot import *
from fullGraph import *
from MCMC import *
from bokeh.models import widgets as wid
from bokeh.plotting import curdoc
from numpy import zeros, reshape, asmatrix, array
import pandas as pd
import bokeh.layouts as layout
from functools import partial
import os

def readPrev(fileName):
    """
    readFile assumes that files are a .txt file with the following format:
     - first an array of the names of all nodes in order followed by a new line character
     - then an array with the most likely edges in it followed by a new lines character
     - then a large matrix with the evolution of all edges followed by a new line character
    """
    f = open(fileName, 'r')
    lines = f.readlines()
    line = lines[0]
    names = line[0:len(line)-2].split(',')
    line = lines[1]
    pOn = line[0:len(line)-2].split(',')
    pOn = list(map(float, pOn))
    line = lines[2]
    pOff = line[0:len(line)-2].split(',')
    pOff = list(map(float, pOff))
    evolve = []
    
    for k in range(3, len(lines)):
        line = lines[k]
        evolve.append(line[0:len(line) -2].split(','))
    
    for i in range(len(evolve)):
        evolve[i] = list(map(float, evolve[i]))
    
    edges = zeros(len(evolve))
    for i in range(len(evolve)):
        cur = evolve[i]
        edges[i] = cur[len(cur) -1]
    size = len(names)
    edges = reshape(asmatrix(edges), (size, size))

    return(evolve, names, edges, pOn, pOff)

def prepCSV(fileName):
    """
    prepCSV allows you to read new data from a csv file, for the purpose of feeding it into 
    the MCMC algorithm
    
    :param fileName: the name of the .csv file to be used
    :return dataDict: a data dictionary with all the necessary components for the algorithm
    """
    data = pd.read_csv(fileName)
    dataDict = {'n': len(data.index), 'C': asarray(data.corr()), 'p': len(data.columns), 'mu': asarray(data.mean()),
                'Sn': asarray(data.cov())}
    names = data.columns
    return(dataDict, names)

def saveMCMC(fileName, names, edgeP):

    """Helps the user by saving the data from their current MCMC Run in case they want to
    pull it up to view again in the future.
    :param fileName: the target destination file to save the new data
    :param names: the names of the nodes
    :param edgeP: the full matrix of edge probabilities to be saved"""
    testSave = open(fileName, 'w+')
    for i in names:
        testSave.write(i + ',') 
        testSave.write('\n')
        exEdge = array(edgeP[0,:])
    for edge in range(len(exEdge[0])):
        toPrint = array(edgeP[:,edge])
        for k in toPrint:
            testSave.write(str(k[0]) + ',')
        testSave.write('\n')
    testSave.close()

def mcmcHelp(hDict):
    
    tab = hDict['tab']; hReq = hDict['selected'].value; opts = hDict['options']
    
    if(hReq == opts[0]):
        tab.child.children[7].children[4] = wid.Div(text = """<p style = 
        "text-align:center;font-family:verdana"> 
        Please selecte a FAQ from the drop down to recieve help on!<p>""")
        
    if(hReq == opts[1]):
         tab.child.children[7].children[4] = wid.Div(text = """<p style = 
         "text-align:center;font-family:verdana"> 
          This tab allows you to either run your own MCMC
          or load the results from a previous MCMC run for 
          analysis.  It defualts to loading an example file
          for you to play around with<p>""")
            
    if(hReq == opts[2]):
         tab.child.children[7].children[4] = wid.Div(text = """<p style = 
         "text-align:center;font-family:verdana"> 
          MCMC is a beautiful sampling algorithm.  When distributions
          are too hard to evaluate directly, either because they are 
          very high dimensional or the integral cannot be evaluated, 
          you have to approximate!  MCMC is a common way to approximate 
          these distribution (for example the one in this application)
          via sampling.  I encourage you to read more about it!<p>""")

def mcmcUpdate(pDict):
    gui = pDict['gui']; fileName = pDict['csv'].value; alpha = float(pDict['alpha'].value)
    burnIn = int(pDict['burnIn'].value); numSteps = int(pDict['numSteps'].value)
    thin = int(pDict['thin'].value); gamma = float(pDict['gamma'].value)
    
    dataDict, names = prepCSV(fileName)
    evolve, gs, lp, _, _ = minIMAP_MCMC(dataDict, burnIn = burnIn, alpha = alpha, justEdges = True, gamma = gamma, 
                 nSteps = numSteps, thinning = thin)
    edges = zeros(len(names)**2)
    for i in range(len(evolve[0,:])):
        cur = array(evolve[:,i]).flatten()
        edges[i] = cur[len(cur)-1]
    size = len(names)
    edges = reshape(asmatrix(edges), (size, size))
    
    gui.tabs[1] = fullTab(edges, names)    
    gui.tabs[2] = edgeTab(edges, names)
    #gui.tabs[3] = learnTab(evolve, names)

prev_val = os.getcwd() + '\\minIMAPapp\\graphExample.txt'
print(prev_val)
csvInput = wid.TextInput(value = "exampleFile2.csv", title = "CSV Data File:")
prevInput = wid.TextInput(value = prev_val, title = "Previous MCMC Session File:")
saveInput = wid.TextInput(value = '', title = "Destination file for saved MCMC data:")

csvButton = wid.Button(label = "Run from CSV", button_type = "success")
prevButton = wid.Button(label = "Show Saved MCMC File", button_type = "success")
saveButton = wid.Button(label = "Save MCMC Info", button_type = "success")


instructions = wid.Div(text = """<p style = "text-align:center;font-family:verdana"> 
    Either run a new MCMC approximation to the data in your CSV
    file by putting in the path to your file and clicking
    run from CSV.  If you simply want to show a previous 
    MCMC run, put the path to the saved MCMC file and click 
    Show Saved MCMC File.  You can also save your Current MCMC 
    run by clicking save after the algorithm has finished. <p>""")
sp0 = wid.Div(width = 60); sp1 = wid.Div(height = 40); 
sp2 = wid.Div(height = 20); sp3 = wid.Div(height = 20) #spacers to make the tab more visually appealing
sp4 = wid.Div(width = 40); sp5 = wid.Div(width = 40)
sp6 = wid.Div(height = 100); sp7 = wid.Div(width = 40)
sp8 = wid.Div(height = 40); sp9 = wid.Div(height = 60)
sp10 = wid.Div(width = 40)

alpha = wid.Slider(start = 0, end = 0.1, value = .05, step = .01, title = "Conditional Independence Cutoff Level")
burnIn = wid.TextInput(value = '20000', title = 'Number of Steps in Burn In:')
numSteps = wid.TextInput(value = '100000', title = 'Number of Steps:')
thinning = wid.TextInput(value = '100', title = 'Steps to wait between samples:')
gamma = wid.TextInput(value = '1', title = "Gamma (penalization in prior):")

helpButton = wid.Button(label = "Help?")
helpOptions = ["", "What should I do?", "What is MCMC?"]
helpDrop = wid.Select(title = "FAQ:", value = helpOptions[0], options = helpOptions)
helpRead = wid.Div(width = 40)
tab0 = wid.Panel(child = layout.row(sp0, layout.column(sp1, instructions), sp10, layout.column(sp4, csvInput, csvButton, sp2, 
                            prevInput, prevButton, sp3, saveInput, saveButton), sp5, layout.column(sp6, alpha, burnIn,
                            numSteps, thinning, gamma), sp7, layout.column(sp8, helpButton, helpDrop, sp9, helpRead)),
                title = "Run MCMC")

helpDict = {'tab': tab0, 'selected': helpDrop, 'options': helpOptions}
helpButton.on_click(partial(mcmcHelp, hDict = helpDict))  

evolve, names, edges, _, _ = readPrev(prev_val)   

tab1 = fullTab(edges, names)    
tab2 = edgeTab(edges, names)
#tab3 = learnTab(evolve, names)
    
imapGUI = wid.Tabs(tabs = [tab0, tab1, tab2])

passDict = {'gui': imapGUI, 'csv': csvInput, "alpha": alpha, 'burnIn': burnIn, 'numSteps': numSteps,
           'thin': thinning, 'gamma': gamma}
csvButton.on_click(partial(mcmcUpdate, pDict = passDict))

curdoc().add_root(imapGUI)
curdoc().title = "Minimial IMAP MCMC"