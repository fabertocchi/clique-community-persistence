from utils import *
import matplotlib.pyplot as plt
import gudhi as gd
import seaborn as sns
import numpy as np

""" toy example to show how the implemented functions work"""

G = nx.Graph()
G.add_edge(0, 1, weight=1)
G.add_edge(1, 2, weight=1)
G.add_edge(2, 3, weight=1)
G.add_edge(0, 3, weight=1)
G.add_edge(1, 3, weight=1)
G.add_edge(0, 4, weight=2)
G.add_edge(3, 4, weight=1)
G.add_edge(4, 5, weight=2)
G.add_edge(4, 6, weight=2)
G.add_edge(5, 6, weight=2)
G.add_edge(0, 5, weight=3)

G_old = G

# to plot and visualize the network
pos = nx.spring_layout(G)
plt.title(f'Graph of toy example network', fontsize=16)
nx.draw(G, pos, with_labels=True)
ws = nx.get_edge_attributes(G, "weight")
nx.draw_networkx_edge_labels(G, pos, edge_labels=ws)
plt.show()

# to extract all cliques from the network G
all_cliques = list(nx.enumerate_all_cliques(G))
# maximal length of cliques
len_max_clique = len(max(all_cliques, key=len))
print("Maximal length of cliques:", len_max_clique)

""" costruction of the k-clique connectivity graphs, considering k = 3 (= len_max_clique)"""
clique_conn_graph, graphs = connectivityGraph(G, len_max_clique)
max_weight = max(w['weight'] for u, v, w in G.edges(data=True))

# to visualize 3-connectivity graphs for all weight from 1 to 3
for w in range(1, len(graphs)):
    g = graphs[w]
    pos = nx.spring_layout(g)
    plt.title(f'3-clique connectivity graph, weight = {w}', fontsize=16)
    nx.draw(g, pos, with_labels=True)
    plt.show()

""" computation of 0-dimensional persistent homology of 3-clique connectivity graph"""

ph = cliquePersistentHomology(graphs, G_old)

# to plot the persistence diagram
print("Persistence diagram: ", ph)
gd.plot_persistence_diagram(persistence=ph)
plt.show()

x = range(0, max_weight+2)
y = []
for i in x:
    y.append(persistenceIndicator(ph, i))

plt.step(x, y, where='post', color="red", linewidth=2)
plt.title(f'Persistence Indicator Function, k = {len_max_clique}')
plt.xlabel('Threshold epsilon')
plt.ylabel('Number of active connected components')
plt.yticks(range(0, 6))
plt.xticks(range(0, max_weight+2))
plt.grid(True, color="grey", linestyle="--", linewidth=0.7)
plt.tight_layout()
plt.show()

""" computation of centrality measure"""

centr = []
for v in G.nodes():
    centr.append((v,cliqueCommunityCentrality(v, G, 5)))

centr_ord = sorted(centr, key=lambda x: (x[1], x[0]), reverse=True)
print("Node with highest centrality measure:", centr_ord[0][0], "value:", centr_ord[0][1])
print("All (nodes, values) ordered from highest centrality measure to lowest:", centr_ord)

