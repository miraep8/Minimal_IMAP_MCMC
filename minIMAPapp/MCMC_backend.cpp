from causaldag import DAG
from causaldag.utils.ci_tests import gauss_ci as gc
from numpy import identity, dot, log, pi, mod, copy, append, asarray, asmatrix, exp, zeros, log1p, concatenate, arange
from numpy.linalg import det
from numpy.random import choice, randint, permutation, uniform

def logLike(DAG, dataDict):
    """
    likelihood essentially return the marginal likelihood of the data given the graph. Based on the 2014 Addendum
    :param DAG: the directed acyclic graph object that we are evaluating
    :param dataDict: the dictionary conatining lots of the relevant information aboout the data, including sample mean
    and covariance.
    
    :return marLik" the marginal likelihood of the current graph structure given the data
    """
    marLik = 0 
    alphaMu = 1 #Currently I arbitrarily fix this hyperparameter - should look more into what this means though
    N = dataDict['n']
    ESS = len(DAG.nodes) + 20 ###  I NEED TO CHANGE THIS - I thought ESS = sample num - num edges estimated
    T = alphaMu*(ESS - len(DAG.nodes) - 1)/(alphaMu + 1)*identity(len(DAG.nodes))  #suggested value for T
    R = T + dataDict['Sn'] + (N*ESS/(N + ESS))*(dot(dataDict['mu'], dataDict['mu']))
    for n in DAG.nodes:
        parents = DAG.upstream(n)
        #p(d|DAG) = prod p(d(Pa_i U {X_i})| DAG)/p(d(Pa_i)| DAG):
        marLik += scoreHelp(parents.union({n}), parents, ESS, R, N, alphaMu, len(DAG.nodes))
    return marLik
    
def scoreHelp(qs, ps, ESS, R, N, aMu, n):
    """
    scoreHelp will compute the log marginal likelihood of a subset of nodes of the graph (typical parents of a
    particular node or the union of the parents of a node with the node itself)
    :param qs: the indices of the set of all parents of node x and x
    :param ps: the indices of the set of all parents of node x
    :param ESS: the effective sample size (also referred to as alpha_w) 
    :param R: matrix which integrate information from the data/prior
    :param N: the number of data points
    :param aMu: a hyperparameter greater than 0
    :param n: the number of nodes in the graph
    
    :return retVal: the marginal likelihood of the data given this part of the graph.
    """
    #define these so that we can properly select the right rows/columns to use:
    xsQ = asarray([[i] for i in qs])
    ysQ = asarray([i for i in qs])
    xsP = asarray([[i] for i in ps])
    ysP = asarray([i for i in ps])
    
    p = len(ps) # number of parents
    
    retVal = (-N/2)*log(pi) + (1/2)*log((aMu)/(N + aMu))
    for i in range(int(N/2)): #rather than using the gamma function, just add the sum of the logs of the remainder
        retVal += log((2*i + ESS - n + p + 1)/2)
    retVal += ((ESS - n + 2*p + 1)/2)*log((aMu*(ESS - n - 1))/(aMu + 1))
    retVal -= ((N + ESS - n + p + 1)/2)*log(det(R[xsQ,ysQ]))
    if not len(ysP) == 0: #determinant of 0x0 matrix is taken to be one
        retVal +=  ((N + ESS - n + p)/2)*log(det(R[xsP,ysP]))
    return(retVal)

def transposePerm(perm):
    """
    returns a random new permutation diven the old permutation, as well as whether there was a full swap
    
    :param perm: the current permutation
    
    :return: A new transpoisiton, two side by side elements are selected at random to be swapped. And min index of transposition 
    """
    newInd = randint(0, len(perm) - 1)
    #Find the neighbor to switch:
    neighbor = mod(newInd + choice([-1,1]), len(perm))
    newPerm = copy(perm)
    newPerm[newInd] = perm[neighbor]
    newPerm[neighbor] = perm[newInd]
    
    return(newPerm, newInd, neighbor)

def newIMAP(perm, alpha, dataDict):
    """
    Creates a new graph object representing the graphical structure of the data, given a permutation.
    :param perm: an array containing an ordering of indices to construct the graph out of
    :param alpha: the confidence level to use for edge construction
    :param dataDict: dictionary structure contianing n: number of data points and C: the correlation matrix of the data
    
    :return newIMAP: a graphical structure representing the new IMAP with the new desired structure. 
    """
    
    IMAP = DAG() #create a new empty graph
    IMAP.add_nodes_from(asarray([k for k in range(len(perm))]))
    
    for i in range(len(perm) -1):
        for j in range(i+1, len(perm)):
            parents = append(asarray(range(i)), asarray(range(i+1,j)))
            S = asarray([perm[p] for p in parents.astype(int)]) #all potential parents of j (other than i) given the perm
            newEdge = gc.gauss_ci_test(dataDict, perm[i], perm[j], S, alpha)
            if newEdge['reject']:
                IMAP.add_arc(perm[i], perm[j])        
    return(IMAP)    
    

def UMI(gCurr, perm, minInd, dataDict, alpha):
    """
    Updates the IMAP in the case that it was not just an exchange of the first and last index (so there are repeated edges).
    
    :param gCurr: DAG structure from the previous acceptance
    :param perm: the new permutation that the new graph is being bult from
    :param minInd: the minimal index exhanged in permutation
    :param dataDict: dictionary structure contianing n: number of data points and C: the correlation matrix of the data
    
    :return gNew: the new graphical structure representing the updated DAG structure.
    """
    
    gNew = DAG()
    edges = gCurr.arcs #take all current arcs from 
    gNew.add_nodes_from(gCurr.nodes)
    if (perm[minInd + 1], perm[minInd]) in edges: #only need to update arcs if the flipped edge was in the original graph
        edges.remove((perm[minInd + 1], perm[minInd]))
        for i in range(minInd):
            for j in asarray([minInd, minInd + 1]):
                parents = append(asarray(range(i)), asarray(range(i+1,j)))
                S = asarray([perm[p] for p in set(parents.astype(int))])
                newEdge = gc.gauss_ci_test(dataDict, perm[i], perm[j], S, alpha)
                if newEdge['reject']: #if we cannot show conditional independnence then we need an edge:
                    edges.add((perm[i], perm[j]))
                elif (perm[i], perm[j]) in edges:
                    edges.remove((perm[i], perm[j]))
        edges.add((perm[minInd], perm[minInd + 1]))
    gNew.add_arcs_from(edges)
    return(gNew)

def minIMAP_MCMC(dataDict, burnIn = 50, alpha = .05, justEdges = False, gamma = 1, nSteps = 500, thinning = 1):
    
    """
    Runs an MCMC algorithm to explore the minimal IMAPs which describe the data.

    :param dataPath: string with the folder address of the data
    :param burnIn: the number of steps to throw out to allow for a burnIn period for the MCMC algorithm. Defaults to 50
    :param alpha: confidence level to use as cutoff for the data.  Defaults to 0.05
    :param justEdges: a boolean value to allow for the user to just choose to get edge estimates rather than all sampled graphs
    :param gamma: the sparsity strength.  Defaults to 1
    :param nsteps: the number of steps to run the MCMC algorithm for.  Defaults to 500
    :param thinning: how many steps to wait before saving a graph.  Defualts to 1

    :return edgeProb: the expectation value of the probability that a given edge is a feature
    :return edgeCount: the evolving counts which show where new edges are seen
    :return logPart: the log partition function
    :return logProb: the matrix which shows the evolution of the log of the posterior for each edge
    :return graphSamples: a collection of all graphs sampled in the MCMC algorithm. 
    """

    perm = permutation(arange(dataDict["p"])) #generate random new permutation to start searching with
    gMap = newIMAP(perm, alpha, dataDict) #create a new random IMAP
    graphSamples = [] # the sampled graphs we return to the user
    count = 0 ; norm = 0 #to help us count number of accepted graphs/number saved graphs
    postOn = [0 for i in range(len(perm)**2)]; postOff = [0 for i in range(len(perm)**2)];
    logProb = []
    edgeCounts = asmatrix([[0 for k in range(len(perm)**2)]])

    for s in range(nSteps):
        if mod(s, int(nSteps/20)) == 0:
            print(str(int(100*s/nSteps)) + '% Done!')
        
        nPerm, i, j = transposePerm(perm)
        if min(i,j) == 0 and max(i,j) == dataDict["p"] - 1: #ie the first and last nodes are exchanged
            gNew = newIMAP(nPerm, alpha, dataDict)
        else:
            gNew = UMI(gMap, nPerm, min(i,j), dataDict, alpha)
        oldLik = logLike(gMap, dataDict)
        newLik = logLike(gNew, dataDict)
        tProb = (newLik - gamma*len(gNew.arcs)) - (oldLik - gamma*len(gMap.arcs))
        if uniform() < exp(tProb): #If np.exp(tProb) > 1 then we will always accept - as required
            count += 1 #by using count rather than s we can gaurantee that we don't skip over the thinning step
            gMap = gNew
            perm = nPerm
            if count > burnIn and mod(count, thinning) == 0:
                norm += 1
                post = newLik - gamma*len(gNew.arcs)
                newProbs = zeros([len(perm)**2])
                for m in range(len(perm)):
                    for n in range(len(perm)):
                        ind = m*len(perm) + n
                        newProbs[ind] = (norm - 1)*edgeCounts[(len(edgeCounts[:,0]) -1), ind]
                        if (m,n) in gMap.arcs:
                            postOn[ind] += post
                            ind = m*len(perm) + n
                            newProbs[ind] += 1
                        else:
                            postOff[ind] += post
                
                edgeCounts = concatenate((edgeCounts, asmatrix(newProbs/norm)), axis = 0)
                if not justEdges:
                    graphSamples.append(gMap)
                logProb.append(post)
    
    return(edgeCounts, graphSamples, logProb, postOn, postOff)