from polygon import Polygon, Point, visualize_polygon, Triangulation, visualize_triangulation
from flip_graph import flip_graph, visualize_graph
import sys, time

def read_input(inputFile):
    """ Read and parse the input file """
    points = []
    with open(inputFile, "r") as f:
        length = int(f.readline())
        for line in f.readlines():
            tokens = line.strip().split()
            new_point = Point(float(tokens[0]), float(tokens[1]))
            points.append(new_point)
    
    return length, points

if __name__ == "__main__":
    filename = sys.argv[1]
    length, points = read_input(filename)
    poly = Polygon(length, points)
    visualize_polygon(poly)
    
    trig = Triangulation()
    output = trig.triangulate(poly)
    for i, o in enumerate(output):
        visualize_triangulation(poly, o, "Polygon Number " + str(i))
    
    # start = time.time()
    # graph = flip_graph(poly)
    # end = time.time()
    # print("Time: ", end - start)
    
    # visualize_graph(graph)