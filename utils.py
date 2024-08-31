import networkx as nx
from networkx.algorithms.community import k_clique_communities
import math
from unionfind import *


def extractEdgesWithWeight(G):
    """ extracts tuples (edge1,edge2,weight) from the given network and returns a list of tuples """
    # NOTICE: it is not needed to sort this tuples since G[u][v]=G[v][u]
    
    elems = []
    for (u,v,w) in G.edges(data=True):
        weight = w["weight"]
        elems.append((u, v, weight))   
     
    return elems

def filtration(G, epsilon):
    """ computes the filtration of given graph, with weight <= epsilon """
    
    G_epsilon = nx.Graph()
    for (u,v,w) in G.edges(data=True):
        weight = w["weight"] # to extract the weight

        # add edges with weight less or equal than epsilon
        if weight <= epsilon:
            G_epsilon.add_edge(u, v, weight=weight)

    return G_epsilon

def cliqueWeight(clique, G):
    """ computes the weight of a single k-clique as the maximum weight of its subsets """
    # NOTICE: to compute the maximum weight of the subset of a 3- or higher dimensional cliques 
    # is equivalent to compute the maximum weight of all the possible edges

    clique = list(clique)
    weight = -math.inf
    
    if len(clique) == 1:    # vertex
        weight = 0

    elif len(clique) == 2:  # edge
        u = clique[0]
        v = clique[1]
        weight = G[u][v]["weight"]

    elif len(clique) >= 3:  # triangle and higher dimensional cliques
        for u in clique:
            for v in clique:
                if u != v:
                    w = G[u][v]["weight"]
                    if weight < w:
                        weight = w
        
    return weight

def weightConnEdges(clique1, clique2, G_k):
    """ computes the weight of edges of a k-clique connectivity graph 
    as the maximum weight between the two k-cliques that are the vertices of that edge """
    
    w1 = cliqueWeight(list(clique1), G_k)
    w2 = cliqueWeight(list(clique2), G_k)

    return max(w1, w2)

def connectivityGraph(G, k, verbose=False):
    """ constructs k-clique connectivity graph of G, with a vertex for every k-clique of G
     and an edge between \sigma and \sigma' if they are adjacent, hence if they share (k-1) vertices """

    graphs = []
    # to find the maximum weight contained in the given network
    max_weight = max(w['weight'] for u, v, w in G.edges(data=True))

    for w in range(max_weight+1):    # iterating over all weights

        if verbose == True:
            print(f"------- w={w} -------")
        clique_conn_graph = nx.Graph()      # initialize an empty graph
        G_w = filtration(G, w)              # compute the w-th filtration

        cliques_w = list(nx.enumerate_all_cliques(G_w))   # to extract all cliques of G_w
        
        for c1 in cliques_w:
            if len(c1) == k: # consider only k-cliques
                c1 = set(c1)
                if tuple(c1) not in clique_conn_graph.nodes():
                    clique_conn_graph.add_node(tuple(c1))   # add a node for every k-clique of G_w
                
                for c2 in cliques_w:
                    if len(c2) == k:    # consider only k-cliques
                        c2 = set(c2)

                        if c1 == c2:
                            continue
                        
                        # to check whether the two considered cliques are adjacent (whether they share k-1 vertices)
                        diff = 0
                        for element in c1:
                            if element not in c2:
                                diff += 1
                        
                        if diff == 1:   # the two cliques share k-1 edges, hence they are adjacent
                            clique_conn_graph.add_edge(tuple(c1), tuple(c2), weight=weightConnEdges(c1, c2, G_w))
                        
        graphs.append(clique_conn_graph)
        if verbose == True:
            print(f"{k}-Clique connectivity graph: edges = {clique_conn_graph.edges()}, nodes = {clique_conn_graph.nodes()}")

    return clique_conn_graph, graphs

def cliquePersistentHomology(graphs, G_old, verbose=False):
    """ computes 0-dimensional persistent homology of k-clique connectivity graph 
    using algorithm 1 (from "Clique Community Persistence: a Topological Visual Analysis Approach for Complex Networks) """
    
    uf = UnionFind()    # initialize empty Union-Find structure
    d = []              # initialize empty persistence diagram
    
    for w in range(len(graphs)):
        if verbose == True:
            print(f"---- w = {w+1} ----")
        elements = extractEdgesWithWeight(graphs[w])
            
        for node in graphs[w].nodes():
            uf.add(node)
    
        if verbose == True:
            print("uf: ", uf)

        for (u, v, weight) in sorted(elements, key=lambda x: (x[2], x[0], x[1])):   # sort edges w.r.t. weight and then vertices
            c1 = uf.find(u)
            c2 = uf.find(v)
            if verbose == True:
                print("weight = ", weight)
                print("u = ", u, "c1 = ", c1)
                print("v = ", v, "c2 = ", c2)

            if c1 == c2:
                continue

            w1 = cliqueWeight(uf[c1], G_old)
            w2 = cliqueWeight(uf[c2], G_old)

            if verbose == True:
                print("uf[c1] = ", uf[c1])
                print("computing w1")
                print(f"w1 = {w1}")
                print("uf[c2] = ", uf[c2])
                print("computing w2")
                print(f"w2 = {w2}")
        
            if w1 <= w2:     # c1 is older component, merge c2 into it
                uf.union(uf[c1], uf[c2])
                if verbose == True:
                    print("uf from if", uf)

                if w2 != weight:
                    d.append((0, (w2, weight)))
                    if verbose == True:
                        print("added", (0, (w2, weight)))
            
            else:           # c2 is older component, merge c1 into it
                uf.union(uf[c2], uf[c1])
                if verbose == True:
                    print("uf from else", uf)

                if w1 != weight:
                    d.append((0, (w1, weight)))
                    if verbose == True:
                        print("added", (0, (w1, weight)))

    if verbose == True:
        print("number of components", len(uf.components()))   

    # to add persistence of components that are never destroyed
    for comp in uf.components():
        weight = []
        for elem in comp:
            weight.append(cliqueWeight(elem, G_old))
        birth = min(weight)
        d.append((0, (birth, math.inf)))
    
    if verbose == True:
        print(f"Persistence diagram: {d}")

    return d

def persistenceIndicator(diagram, epsilon):
    """ computes number of active connected components for given threshold epsilon """

    n = 0
    for comp in diagram:
        birth = comp[1][0]
        death = comp[1][1]
        if epsilon >= birth and epsilon < death:
            n += 1

    return n

def glyph(diag, divider, length):
    """ condensed glyph that gives a summary of clique community activity """
    
    steps = length // divider   # discretize the domain into uniformly spaced bins
    values = []
    for step in range(steps):
        value = -math.inf
        for epsilon in range(step * divider, step * divider + divider+1):
            value = max(value, persistenceIndicator(diag, epsilon)) # maximum value of persistence indicator in that range
        values.append(value)
        
    return values

def computePersistence(graphs, G_old):
    """ keeps track of each community and its persistence """
    
    uf = UnionFind()
    d_withcomp = []
    
    for w in range(len(graphs)):
        elements = extractEdgesWithWeight(graphs[w])
    
        for node in graphs[w].nodes():
            uf.add(node)

        for (u, v, weight) in sorted(elements, key=lambda x: (x[2], x[0], x[1])):
            c1 = uf.find(u)
            c2 = uf.find(v)

            if c1 == c2:
                continue

            w1 = cliqueWeight(uf[c1], G_old)
            w2 = cliqueWeight(uf[c2], G_old)

            
            if w1 < w2:
                uf.union(uf[c1], uf[c2])

                if w2 != weight:
                    d_withcomp.append((uf.component(v), (w2, weight)))    # component and its persistence
            
            else:
                uf.union(uf[c2], uf[c1])

                if w1 != weight:
                    d_withcomp.append((uf.component(u), (w1, weight)))    # component and its persistence
  
    # components that are never destroyed
    for comp in uf.components():
        weight = []
        for elem in comp:
            weight.append(cliqueWeight(elem, G_old))
        birth = min(weight)
        d_withcomp.append((set(comp), (birth, math.inf)))
    
    return d_withcomp

def cliqueCommunityCentrality(v, G, weight_infinity=32):
    """ measures the relevance of vertex v considering persistence of all communities v is part of """
    all_cliques = list(nx.enumerate_all_cliques(G))
    len_max_clique = len(max(all_cliques, key=len))
    centrality = 0
    ph = []
    
    # to compute persistence of all communities for all k 
    for k in range(2, len_max_clique+1):
        clique_conn_graph, graphs = connectivityGraph(G, k)
        ph.append(computePersistence(graphs, G))

    for element in ph:
        for el in element:
            for e in el[0]:
                if v in e:  # check whether v is in the considered community
                    birth = el[1][0]
                    death = el[1][1]

                    #  extract persistence of the community
                    persistence = abs(death-birth) if death != float('inf') else weight_infinity
                    centrality += persistence
                    break

    return centrality
