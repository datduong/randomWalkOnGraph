import igraph as ig 
import numpy as np 
import scipy.stats as stat 
import matplotlib.pyplot as plt
# import cairo as cairo

class randNetwk (ig.Graph): # random network inherits from Graph 
    def __init__ (self):
        ig.Graph.__init__(self) 
        self.network = None ## keep the actual network. 

    def keep_edges (self,prob): 
        """ make an edge with exact prob Prob P(x=1)=Prob, x=1 implies connection """
        """ keep an edge with some prob """
        self.prob = prob 
        if (self.network is None): 
            print "must create edges"
        else: 
            lenEdges = len(self.network.get_edgelist()) # number of edges 
            randprob = np.random.binomial(1,prob,lenEdges) # draw a 1 with prob=p. 
            edgesIdx = np.linspace(0,lenEdges-1,num=lenEdges) # index of edges 
            edgesIdx = edgesIdx [ randprob == 0 ] # edges to drop 
            edgesIdx = map(int,edgesIdx.tolist()) # convert to int. 
            self.network.delete_edges(edgesIdx) # drop edges 
           
    def make_full_edges (self, n): 
        self.network = ig.Graph.Full(n); # n is number of nodes
                        
    def return_degree  (self): 
        self.hist = self.network.degree()
        
    def hist_degree (self,kind) :
        
        self.return_degree()
        font = {'family' : 'serif', 'color'  : 'darkred', 'weight' : 'normal', 'size' : 14 }
        
        if kind == 'prob':
            plt.hist(self.hist) 
            plt.title('Degree distribution p='+str(self.prob), fontdict=font)
        if kind == 'grow':
            plt.hist(self.hist) 
            plt.title('Degree distribution of a growing network ', fontdict=font)
        if kind == 'in' :
            self.hist = self.network.indegree()  
            plt.hist(self.hist) 
            plt.title('in_Degree distribution', fontdict=font)
        if kind == 'out': 
            self.hist = self.network.outdegree() 
            plt.hist(self.hist) 
            plt.title('out_Degree distribution', fontdict=font)
            
        plt.xlabel('degree', fontdict=font)
        plt.ylabel('count', fontdict=font)
        
    def make_plot (self,lay_out): # see reference for lay_out 
        vis = {} # visual input 
        vis["vertex_size"]=1
        vis["edge_width"]=.1
        self.netwkPlot = ig.plot(self.network, layout = lay_out, **vis) ## pass in vis. 
        # to see plot, do: rand.netwkPlot.show() 
   
    def greedy_subNetwk (self):
        self.greedyNetwk = rand.network.community_fastgreedy() # returns dendrogram 
        self.greedyNetwkClst = self.greedyNetwk.as_clustering() # break into optimal clusters 
        print 'optimal number of clusters for this dendrogram ' + str(self.greedyNetwk.optimal_count)
        print "the modularity score for community network" + str(self.greedyNetwkClst.modularity)
                                                                                     
    def connected (self): 
        print 'connected: ' + str(self.network.is_connected())    
       
    def begin_growing(self,N): 
        """ add a node to existing graph, and make a connection of this new node to old ones with prob f(new,old)"""
        # N: number of nodes at end of run-time
        n = 1 # start n
        self.network = ig.Graph(1) # add one node to begin 
        while (n<=N): 
            sumDegree = np.sum(self.network.degree()) # current degree 
            prob = [] # start prob     
            for i in range(0,n):       
                if sumDegree == 0: 
                    prob.append(1.0/(n+1)) ## only happens if sumDegree=0, n=1 then 1/2 connect 1/2 not, n=2 then 1/3 connect to 1, 1/3 connect to 2, or 1/3 not 
                else:
                    if ( self.network.degree()[i] == 0 ): ## possible that one old node has no connection at all (due to being unlucky)
                        prob.append( 0.05 ) ## ensure will be connected at some point 
                    else: 
                        prob.append( 1.0 * self.network.degree()[i] /sumDegree )   
            ''' add a node '''
            self.network.add_vertex() 
            prob = np.array(prob)     
            ''' add edges '''
            makeConnect = stat.bernoulli.rvs(prob) # if 1, then connect 
            if isinstance(makeConnect,int):
                makeConnect = np.array ( [makeConnect] ) 
            for i in range(0,n):
                if makeConnect[i] == 1:
                    self.network.add_edge(source=n,target=i)    
            # at start, n=1, , so now increase number of n 
            n += 1
			
    def forest_fire (self, param):
        try: 
            self.network = ig.Graph.Forest_Fire(**param) ## pass a bunch of stuffs in as a dictionary 
            print 'diameter' + str (self.network.get_diameter())
        except: 
            print 'parameters not given'
            
    def giant_connected (self): 
        if self.network is None: 
            print 'must make a network' 
        self.giantcc = self.network.clusters()
        self.giantcc = self.giantcc.giant() # giant component may not be unique 
        print 'nodes in giant connected component: ' + str(len(self.giantcc.degree()))
        self.ratioGC_trueGraph = len(self.giantcc.degree())*1.0 / len(self.network.degree())
        print 'ratio of |Gcc| over |G| ' + str(self.ratioGC_trueGraph)
        
    def randomWalk (self,stop): 
        if (self.network is None): 
            print "must create the network"     
        ''' start random position, walk until @stop '''
        startNode = np.random.random_integers(0, len(rand.network.degree())-1 ) # start position 
        t = 1
        nodeTravelto = np.array([startNode]) # node traveled to 
        distTravel = np.array([0]) # distance travel 
        while (t <= stop ):
            ''' randomly chose a neighbor '''
            neighbors = self.network.neighbors(startNode) # get the neighbor nodes. 
            if (neighbors.__len__()==0): # graph is disconnected 
                neighbors = self.network.neighborhood(startNode)
    
            # self.network.degree(neighbors)
            currentNode = np.random.choice(neighbors) # randomly choose a neighbor 
            nodeTravelto = np.append(nodeTravelto,currentNode)
            ''' find shortest path from current node to node zero '''
            dist = self.network.shortest_paths(nodeTravelto[0],currentNode) 
            distTravel = np.append (distTravel,dist)
            ''' update info for next round '''
            startNode = currentNode # update location 
            t = t+1 # could use for-loop instead. whatever.
            
        ''' return both arrays '''
        return nodeTravelto, distTravel
            
    def manyRandomWalks (self,numerate,stop):
        if (self.network is None): 
            print "must create the network"     
        ''' @numerate is the number of walkers -- the walking dead '''
        self.randomWalkersNode = np.ones((numerate,stop+1)) # plus one due to stop=10, then there are 11 total (due to starting location)
        self.randomWalkersDist = np.ones((numerate,stop+1))
        for i in range(0,numerate): # run for many walkers 
            nodeTravelto, distTravel = self.randomWalk(stop)
            self.randomWalkersNode[i,] = nodeTravelto
            self.randomWalkersDist[i,] = distTravel
        ''' the mean & var distance walk '''    
        self.meanDistWalk = np.mean(self.randomWalkersDist,0)
        self.varDistWalk = np.var( np.diff ( self.randomWalkersDist ) , 0)
        self.varDistWalk = self.varDistWalk.cumsum()
        
    def powerLawGraph (self,n,m, degree): # m number of edges in the graph, n is num-node 
        # self.network = ig.Graph.Barabasi(n, m, outpref=False, directed=False, power=3, zero_appeal=1, implementation="psumtree", start_from=None)
        self.network = ig.Graph.Static_Power_Law(n ,m ,exponent_out = degree);
		# returns a directed or undirected graph with the prescribed power-law degree distributions.

    def plotMeanDist (self):  
        p1=plt.plot(self.meanDistWalk)
        font = {'family' : 'serif', 'color'  : 'darkred', 'weight' : 'normal', 'size' : 14 }
        plt.title('mean distance traveled', fontdict=font) 
        plt.xlabel('time', fontdict=font)
        plt.ylabel('average shortest distance', fontdict=font)             
        plt.legend([p1[0]], ['Mean'], 'best', numpoints=1)
        print ('diameter ' + str(self.network.diameter()) )  

    def plotVarDist (self):  
        p2=plt.plot(self.varDistWalk)
        font = {'family' : 'serif', 'color'  : 'darkred', 'weight' : 'normal', 'size' : 14 }
        plt.title('variance distance traveled', fontdict=font) 
        plt.xlabel('time', fontdict=font)
        plt.ylabel('distance', fontdict=font)             
        plt.legend([p2[0]], ['Var'], 'best', numpoints=1) 

    def plotHistFinalNode (self):    
        font = {'family' : 'serif', 'color'  : 'darkred', 'weight' : 'normal', 'size' : 14 }
        plt.title('histogram of degree at final node', fontdict=font) 
        plt.xlabel('degree', fontdict=font)
        plt.ylabel('count', fontdict=font)             
        plt.hist(rand.network.degree( map( int , self.randomWalkersNode[:,stop]) ) )
 
 
##############################################################################
           
'''randomly connect nodes with prob=p '''

rand = randNetwk()                                                                     
rand.make_full_edges(1000) ## create a undirected, complete graph with 1000 nodes 
rand.keep_edges(.01) ## keeping edge with some probability, else delete the edge  
rand.connected()       

stop=100 # stop time. 
rand.manyRandomWalks(100,stop) ## 100 walkers, each starts at a random node in the graph

rand.plotMeanDist() ## the mean distance from starting point of the 100 random walks
rand.plotVarDist() ## the variance of distance from starting point of the 100 random walks
rand.plotHistFinalNode() ## the degree of the final node (destination) of the 100 random walks 

rand.make_plot('kk') ## to see plot, 
rand.netwkPlot.show() ## time consuming 
# rand.greedy_subNetwk()
# rand.giant_connected()

##############################################################################

''' make barabasi graph, then do random walk '''  

rand = randNetwk()  
rand.powerLawGraph(1000,999,3) ## 1000 nodes, 999 edges. In this undirected graph, new node i, connects to existing j with p_{ij} propto to degree(j)/ sum_k degree(k). 
## degree distribution is ~k^{-3}, where 3 is an input to the function.     

stop = 200
rand.manyRandomWalks(100,stop) ## 100 walkers
rand.plotMeanDist()
rand.plotVarDist()
rand.plotHistFinalNode()
