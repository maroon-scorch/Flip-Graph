from polygon import *
import networkx as nx
from scipy import sparse
import matplotlib.pyplot as plt
import numpy as np
import time
from flip_graph import flip_graph, visualize_graph

# https://www.geeksforgeeks.org/power-method-determine-largest-eigenvalue-and-eigenvector-in-python/
def power_iteration(A, max_iter: int):
    # Ideally choose a random vector
    # To decrease the chance that our vector
    # Is orthogonal to the eigenvector
    x = np.random.rand(A.shape[1])
    tol = 1e-6
    lam_prev = 0
    
    for i in range(max_iter):
        # Compute the updated approximation for the eigenvector
        x = A @ x / np.linalg.norm(A @ x)
    
        # Compute the updated approximation for the largest eigenvalue
        lam = (x.T @ A @ x) / (x.T @ x)
    
        # Check if the approximations have converged
        if np.abs(lam - lam_prev) < tol:
            break
    
        # Store the current approximation for the largest eigenvalue
        lam_prev = lam

    return lam

if __name__ == "__main__":
    
    start = time.time()

    iter = 10
    iter2 = -1
    x = np.arange(4, iter2 + 1)
    # x = np.arange(5, iter2 + 1);
    
    y = [];
    
    for i in range(4, iter + 1):
        print("Polygon ", i);
        polygon = regular_polygon(i)
        graph = flip_graph(polygon)
        
        G = nx.Graph(graph)
        result = nx.algebraic_connectivity(G)
        # L = nx.laplacian_matrix(G).toarray()
        
        # reduced_L = np.delete(L, np.s_[-1:], axis=1)
        # reduced_L = np.delete(reduced_L, np.s_[-1:], axis=0)
        
        # val = power_iteration(L, 200)
        # eigenvalues, eigenvectors = np.linalg.eigh(L)

        print(result)
        # print(G.number_of_edges())
        # print(eigenvalues[1])
        # print(eigenvalues[2])
        # print(eigenvalues[3])
        # # print(nx.diameter(G))
        y.append(result);
        print("-------------------------------------")
        
    for i in range(iter + 1, iter2 + 1):
        print("Polygon ", i);
        polygon = regular_polygon(i)
        graph = flip_graph(polygon)
        
        G = nx.Graph(graph)
        degree_dict = nx.degree(G)
        degree_list = [x[1] for x in degree_dict]
        lap_matrix = sparse.diags(degree_list, 0)-nx.adjacency_matrix(G)
        eigval, eigvec = sparse.linalg.eigsh(lap_matrix, 3, sigma=0, which='LM')

        print(G.number_of_edges())
        print(eigval)
        print(eigval[1])
        y.append(eigval[1]);
        # print(nx.diameter(G))
        print("-------------------------------------")
    
    y = np.asarray(y);
    plt.plot(x, y, 'bo--', linewidth=2, markersize=12)
    # plt.show()
    plt.savefig('plot.png')
    
    end = time.time()
    print("Time: ", end - start)