from polygon import *
import networkx as nx
import matplotlib.pyplot as plt
import copy, itertools


trig = Triangulation()

def visualize_graph(graph):
    """ Given a graph in the vertex dictionary representation, graphs it """
    G = nx.Graph(graph)
    # Displaying the graph
    pos = nx.spring_layout(G)
    # Drawing the graph
    nx.draw_networkx_nodes(G, pos, node_size=700)
    nx.draw_networkx_edges(G, pos, width=6)
    
    nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")
    
    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    plt.tight_layout()
    plt.show()


def dfs(graph, visited, n, path, triangles):
    current_vertex = path[-1]
    visited[current_vertex] = True
    
    # We have terminated our steps
    if n == 0:
        visited[current_vertex] = False
        # Check the current vertex can reach start
        if path[0] in graph[current_vertex]:
            triangles.append(path)
            
    # Search pathes of length n - 1 from current vertex
    for i in range(len(graph.keys())):
        if visited[i] == False and i in graph[current_vertex]:
            new_path = path.copy() + [i]
            dfs(graph, visited, n - 1, new_path, triangles)
    
    visited[current_vertex] = False
    
    
def diagonals_to_triangles(diagonals, n):
    """ Given a list of diagonals on a polygon that forms the triangulation,
    returns the list of triangles divided by this. """
    
    # print(diagonals)
    # Creates a graph based on the information
    graph = {}
    # Maybe we could switch set to list later (no uniqueness issues)
    for i in range(n):
        graph[i] = set([(i - 1) % n, (i + 1) % n])
    for i, j in diagonals:
        graph[i].add(j)
        graph[j].add(i)
    # print(graph)
    
    # A triangle is equivalent to a cycle of length 3 in this graph
    triangles = []
    visited = [False] * n
    for i in range(n - 2):
        dfs(graph, visited, 2, [i], triangles)
        # Every triangle with index i has been enumerated already
        visited[i] = True
    
    for t in triangles:
        t.sort()
        
    assert len(triangles) % 2 == 0
    
    return triangles[::2], graph
    

def find_quadrilateral(d, triangles, n, graph):
    """ Finds the quadrilateral that contains this diagonal d """
    # Since this is an diagonal of the quadrilateral, the index of the other two points must be adjacent
    # to the two points of the diagonal
    
    index = graph[d[0]].union(graph[d[1]])
    
    for i, j in itertools.combinations(index, 2):
        t1 = [d[0], d[1], i]
        t2 = [d[0], d[1], j]
        t1.sort()
        t2.sort()
        # print(t1)
        # print(t2)
        # print("------------------")
        
        if t1 in triangles and t2 in triangles:
            return [i, j, d[0], d[1]]
    print(d)
    print(triangles)
    print("Shouldn't reach here")
    

def has_edge(trig_1, trig_2, n, triangles, graph):
    """ Determines if there's an edge between Trig 1 and Trig 2"""

    # print(trig_1)
    # print(trig_2)
    
    for d in trig_1:
        i, j, _, _ = find_quadrilateral(d, triangles, n, graph)
        t1 = copy.deepcopy(trig_1)
        t2 = copy.deepcopy(trig_2)
        t1.remove(d)
        
        if i < j:
            t1.append((i, j))
        else:
            t1.append((j, i))
        
        t1.sort(key = lambda x: x[0] + (1/n)*x[1])
        if t1 == t2:
            return True
    
    return False

def flip_graph(poly):
    """ Generates the flip graph of a given polygon """
    # Find all triangulates of the polygon:
    graph = {}
    
    # These diagonals are specified by points
    all_triangulations = trig.triangulate(poly)
    length = len(all_triangulations)
    
    # These diagonals are specified by index
    index_triangulations  = []
    for tri in all_triangulations:
        input = map(lambda ed: (poly.vertice.index(ed[0]), poly.vertice.index(ed[1])), tri)
        input = list(map(lambda ed: (ed[1], ed[0]) if ed[0] > ed[1] else ed, input))
        input.sort(key = lambda x: x[0] + (1/poly.n)*x[1])
        index_triangulations.append(input)
    
    print(index_triangulations)
    print("Finsihed Triangulation")
    print("Number of Triangulations: ", len(index_triangulations))
        
    for i in range(length):
        graph[i] = []
        
    for i in range(length):
        print(i)
        triangles, g_dict = diagonals_to_triangles(index_triangulations[i], poly.n)
        for j in range(i + 1, length):
            if has_edge(index_triangulations[i], index_triangulations[j], poly.n, triangles, g_dict):
                graph[i].append(j)
                graph[j].append(i)
    
    print("Finished Making Graph")
    
    return graph

if __name__ == "__main__":
    # polygon = regular_polygon(5)
    # polygon = Polygon(6, [Point(0, 0), Point(1, 0), Point(2, 0), Point(2, 1), Point(1, 1), Point(0, 1)])
    polygon = Polygon(5, [Point(0, 0), Point(2, 0), Point(2, 2), Point(1, 1), Point(0, 2)])
    visualize_polygon(polygon)
    graph = flip_graph(polygon)
    print(graph)
    visualize_graph(graph)
    
    # ----------------------------------
    # poly = regular_polygon(6)
    # all_triangulations = trig.triangulate(poly)
    # print(all_triangulations[0])
    # index_triangulations  = []
    # for tri in all_triangulations:
    #     input = map(lambda ed: (poly.vertice.index(ed[0]), poly.vertice.index(ed[1])), tri)
    #     input = list(map(lambda ed: (ed[1], ed[0]) if ed[0] > ed[1] else ed, input))
    #     input.sort(key = lambda x: x[0] + (1/poly.n)*x[1])
    #     index_triangulations.append(input)
    # print(has_edge(index_triangulations[0], index_triangulations[1], poly))
     
        
    # diagonals_to_triangles(index_triangulations[0], 8)
    # visualize_triangulation(polygon, all_triangulations[0])