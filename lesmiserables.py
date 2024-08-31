from utils import *
import matplotlib.pyplot as plt
import gudhi as gd
import seaborn as sns
import numpy as np

""" all functions applied to Les Miserables dataset """

# to import Les Miserables dataset and construct the corresponding weighted network
file_path = './lesmiserables_num.gml'
G = nx.read_gml(file_path)
G_old = G.copy()

# to compute maximal weight contained in the network
max_weight = max(w['weight'] for u, v, w in G.edges(data=True))

# to invert the weights
# NOTICE: edge weights in dataset correspond to number of co-occurrences between two characters
# => invert weights since we consider edge weights to correspond to proximity
for (u,v,w) in G_old.edges(data=True):
    weight = max_weight - w["weight"]
    G[u][v]["weight"] = weight

# new maximal weight after inversion
max_weight = max(w['weight'] for u, v, w in G.edges(data=True))

# to visualize the network 
# pos = nx.spring_layout(G)
# nx.draw(G, pos, with_labels=True)
# ws = nx.get_edge_attributes(G, "weight")
# nx.draw_networkx_edge_labels(G, pos, edge_labels=ws)
# plt.show()


# to understand better the structure of the network
print(len(G.edges(data=True)), "edges")  # 254 edges
print(len(G.nodes()), "nodes")   # 77 nodes

# to extract all cliques of G
all_cliques = list(nx.enumerate_all_cliques(G))

# length of clique with maximal length
len_max_clique = len(max(all_cliques, key=len))
print("Length of clique with maximal length:", len_max_clique)

""" Considering the special case of k = 4 """

# compute the 4-clique connectivity graph
clique_conn_graph, graphs = connectivityGraph(G, 4)

# to plot the evolution of 4-clique connectivity graph w.r.t. weight threshold
for w, graph in enumerate(graphs):
    if w > 15:
        plt.figure()
        pos = nx.spring_layout(graph)
        plt.title(f"Weight = {w}", fontsize=16)
        nx.draw(graph, pos, with_labels=True)
        plt.show()

# to check number of communities (for weight=30)
communities = list(k_clique_communities(G, 4))
print("Communities:", communities)
print("Number of communities:", len(communities))
print("Connected components from graph:", nx.number_connected_components(clique_conn_graph))


""" Considering all k-clique connectivity graphs"""

# construct k-clique connectivity graphs and compute their 0-dimesional persisten homology
ph = []
for k in range(2, len_max_clique+1):
    clique_conn_graph, graphs = connectivityGraph(G, k)

    pos = nx.spring_layout(clique_conn_graph)
    plt.title(f"{k}-clique connectivity graph", fontsize=16)
    nx.draw(clique_conn_graph, pos, with_labels=True)
    plt.show()

    # to check the number of communities
    communities = list(k_clique_communities(G, k))
    print(f"k = {k}, number of communities: {len(communities)}")
    print("Connected components from graph:", nx.number_connected_components(clique_conn_graph))

    ph_k = cliquePersistentHomology(graphs, G)

    # create a list with all persistence pairs for all k to construct the persistence diagram
    for el in ph_k:
        ph.append(el)

# to plot the persistence diagram
print("persistence diagram:", ph)
gd.plot_persistence_diagram(persistence=ph)
plt.show()

""" Persistence indicator function and clique community centrality """

# compute persistence diagrams for every k
ph_all = []
for k in range(2, len_max_clique+1):
    clique_conn_graph, graphs = connectivityGraph(G, k)
    ph_k = cliquePersistentHomology(graphs, G)
    ph_all.append(ph_k)
    print("k =", k, ph_k)
print(ph_all)

# plot the persistence indicator function (is a step function) for each k separately
x = range(0, max_weight+2)
for k, diag in enumerate(ph_all):
    y = []
    for i in x:
        y.append(persistenceIndicator(diag, i))

    plt.step(x, y, where='post', color="red", linewidth=2)
    plt.title(f'Persistence Indicator Function, k = {k+2}')
    plt.xlabel('Threshold epsilon')
    plt.ylabel('Number of active connected components')
    plt.yticks(range(0, 11))
    plt.xticks(range(0, max_weight+2))
    plt.grid(True, color="grey", linestyle="--", linewidth=0.7)
    plt.tight_layout()
    plt.show()

# # to plot all graphs in one
fig, axes = plt.subplots(nrows=len(ph_all), ncols=1)

for k, diag in enumerate(ph_all):
    y = []
    for xx in x:
        y.append(persistenceIndicator(diag, xx))
    axes[len(ph_all)-k-1].step(x, y, where='post', color="red", linewidth=2)
    axes[len(ph_all)-k-1].set_xticks([])
    axes[len(ph_all)-k-1].set_yticks([])
    if k == 0:
        axes[len(ph_all)-k-1].set_xticks(range(0, 32))

fig.suptitle('Persistence Indicator function', fontsize=16)
plt.tight_layout()
plt.show()

# to visualize the glyph
data = np.zeros((len(ph_all), 8))
for i in range(len(ph_all)):
    data[len(ph_all)-i-1] = glyph(ph_all[i], 4, 32)
print(data)

plt.figure(figsize=(6, 6))
ax = sns.heatmap(data, annot=False, cmap='Blues', cbar=False, square=True, linewidths=0.5)
ax.set(xlabel="Weights ranges", ylabel="k")
ax.xaxis.set_ticklabels([])
ax.xaxis.set_ticks([])
ax.yaxis.set_ticklabels(reversed(range(2, 11)))
plt.show()

# compute the centrality measure for all vertices of G, using clique community centrality measure
centrality = []
for v in G.nodes():
    centrality.append((v, cliqueCommunityCentrality(v, G)))

# sort all pairs (node, value) w.r.t. value in descending order
centrality_ord = sorted(centrality, key=lambda x: (x[1], x[0]), reverse=True)
print("Node with highest centrality measure:", centrality_ord[0][0], "value:", centrality_ord[0][1])
print("All (nodes, values) ordered from highest centrality measure to lowest:", centrality_ord)
print("First five nodes with highest centrality measure:", centrality_ord[:5] )
