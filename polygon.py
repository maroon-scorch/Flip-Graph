import numpy as np
import matplotlib.pyplot as plt
import math
from shapely import geometry
import time

# Class representing a point in R^2
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.vec = [x, y]
        self.array = np.asarray([x, y])

    def __repr__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"

    def __eq__(self, other):
        if type(other) != Point:
            return False
        # return approx_equal(self, other)
        return self.x == other.x and self.y == other.y

    def __hash__(self):
      return hash((self.x, self.y))

def to_point(lst):
    """ Converts a given list to a point """
    assert len(lst) == 2
    return Point(lst[0], lst[1])

def dist(p1, p2):
    """ Given 2 points, find the distance between them """
    return np.linalg.norm(p2.array - p1.array)

def midpoint(p1, p2):
    return to_point((p1.array + p2.array)/2)

# CLass representing a polygon with n vertices
class Polygon:
    def __init__(self, n, vertice):
        # Checks that polygon is not degenerate
        assert len(set(vertice)) == n
        assert n > 2
        # TODO: Check polygon does not have colinear 3 adjacent points
        
        # The number of vertices
        self.n = n 
        # The points of the vertices, given in order of traversal
        self.vertice = vertice
        
        # The edges of the polygon, still in order of traversal
        edges = []
        for idx in range(len(vertice)):
            if idx != len(vertice) - 1:
                edges.append((vertice[idx], vertice[idx + 1]))
            else:
                edges.append((vertice[idx], vertice[0]))
        self.edges = edges
        
        input = list(map(lambda pt: pt.vec, vertice))
        self.shapely_polygon = geometry.Polygon(input)

    def __repr__(self):
        return str(self.vertice)

    def __eq__(self, other):
        if type(other) != Polygon:
            return False
        # If the polygons are simple, and the inputs fulfills the type promise
        # This should be enough (TODO: Check this)
        return self.n == other.n and set(self.vertice) == set(other.vertice)

    def __hash__(self):
      return hash(self.shapely_polygon)

def int_between(x, y):
    if x < y:
        return range(math.ceil(x), math.floor(y) + 1)
    else:
        return range(math.ceil(y), math.floor(x) + 1)
 
def visualize_polygon(poly):
    """ Visualizes a given polygon """
    input = list(map(lambda pt: pt.vec, poly.vertice))
    input.append(input[0])
    
    x_pts, y_pts = zip(*input)
    fig = plt.figure()
    plt.plot(x_pts, y_pts)
    
    for i in range(0, len(x_pts) - 1):
        plt.plot(x_pts[i], y_pts[i], 'r-o')
    
    ax = fig.gca()
    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim() 
    ax.set_xticks(int_between(xmin, xmax))
    ax.set_yticks(int_between(ymin, ymax))
    plt.grid()
    # Title
    plt.title("Polygon")

    plt.show()

def is_line_segment_in_polygon(segment, poly):
    line = geometry.LineString((segment[0].vec, segment[1].vec))
    polygon = poly.shapely_polygon
    boundary = geometry.LineString(list(polygon.exterior.coords))
    intersections = boundary.intersection(line)
    
    # print(intersections)
    # The line can lie on the boundary    
    if intersections.union(line) == intersections:
        return True
    
    # The end points of the line can lie on boundary, but the interior has to be
    # in the polygon
    temp = intersections.difference(geometry.Point(segment[0].x, segment[0].y))
    temp = temp.difference(geometry.Point(segment[1].x, segment[1].y))
    
    if not temp.is_empty:
        # This means that the line definitely crossed the boundary
        return False 
    else:
        # Check if mid point of the line is contained in the Polygon
        mid = midpoint(segment[0], segment[1])
        shapely_mid = geometry.Point(mid.x, mid.y)

        # If it is, then it can't intersect on the boundary because temp is empty
        if not polygon.contains(shapely_mid):
            # If this midpoint is outside the polygon, this is false
            return False
        else:
            # The midpoint is inside the polygon, we can check if both end points are inside
            return (boundary.contains(geometry.Point(segment[0].x, segment[0].y)) and boundary.contains(geometry.Point(segment[1].x, segment[1].y))) or (polygon.contains(geometry.Point(segment[0].x, segment[0].y)) and polygon.contains(geometry.Point(segment[1].x, segment[1].y)))
    
    
    # Else the line can still have 
    # intersections = intersections.difference(line)
    # # print(intersections)
    # if intersections.is_empty:
    #     # Want both end points to be contained in the boundary
    #     return (boundary.contains(geometry.Point(segment[0].x, segment[0].y)) and boundary.contains(geometry.Point(segment[1].x, segment[1].y))) or (polygon.contains(geometry.Point(segment[0].x, segment[0].y)) and polygon.contains(geometry.Point(segment[1].x, segment[1].y)))
    # else:
    #     # This means a point other than the end points lies on the boundary
    #     return False

def find_all_triangulations(poly):
    """ Given a polygon, returns the list of all possible triangulations of the polygon,
    each triangulation is represented by a list of diagonals."""
    # This is already a triangle
    if poly.n == 3:
        return [[]]
    
    output = []
    # This the edge from index 0 to index 1
    fixed_base = poly.edges[0]
    
    for i in range(2, poly.n):
        current_vertex = poly.vertice[i]
        line_1 = (fixed_base[0], current_vertex)
        line_2 = (fixed_base[1], current_vertex)
        
        if is_line_segment_in_polygon(line_1, poly) and is_line_segment_in_polygon(line_2, poly):
            # This means that this is a valid triangle!
            diagonals = [line_1, line_2]
            # Create the two polygons divided by this triangle
            left_index = list(range(1, i + 1))
            right_index = list(range(i, poly.n)) + [0]
            
            left_triangulations = []
            right_triangulations = []
            
            # This means a left polygon exists
            if len(left_index) > 2:
                left_vertices = list(map(lambda i: poly.vertice[i], left_index))
                left_polygon = Polygon(len(left_index), left_vertices)
                left_triangulations = find_all_triangulations(left_polygon)
            # This means a right polygon exists
            if len(right_index) > 2:
                right_vertices = list(map(lambda i: poly.vertice[i], right_index))
                right_polygon = Polygon(len(right_index), right_vertices)
                right_triangulations = find_all_triangulations(right_polygon)
            
            if left_triangulations != [] and right_triangulations != []:
                # Combining the Triangulations
                for left_diagonals in left_triangulations:
                    for right_diagonals in right_triangulations:
                        possible_triangulation = left_diagonals + diagonals + right_diagonals
                        output.append(possible_triangulation)
            elif left_triangulations != [] and right_triangulations == []:
                for left_diagonals in left_triangulations:
                    possible_triangulation = left_diagonals + diagonals
                    output.append(possible_triangulation)
            elif right_triangulations != [] and left_triangulations == []:
                for right_diagonals in right_triangulations:
                    possible_triangulation = diagonals + right_diagonals
                    output.append(possible_triangulation)
            else:
                # Shouldn't reach this case
                print("Shouldn't reach this case")
                    
    return output

def visualize_triangulation(poly, diagonals, title):
    input = list(map(lambda pt: pt.vec, poly.vertice))
    input.append(input[0])
    
    x_pts, y_pts = zip(*input)
    fig = plt.figure()
    plt.plot(x_pts, y_pts)
    
    for i in range(0, len(x_pts) - 1):
        plt.plot(x_pts[i], y_pts[i], 'r-o')
    
    for i, ed in enumerate(diagonals):
        start = ed[0]
        end = ed[1]
        plt.plot([start.x, end.x], [start.y, end.y], 'k-')    
    
    ax = fig.gca()
    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim() 
    ax.set_xticks(int_between(xmin, xmax))
    ax.set_yticks(int_between(ymin, ymax))
    plt.grid()
    # Title
    plt.title(title)

    plt.show()
    
def regular_polygon(n):
    """ Creates a regular polygon """
    vertice = []
    for i in range(n):
        x = math.cos(2*math.pi*i/n)
        y = math.sin(2*math.pi*i/n)
        vertice.append(Point(x, y))
    return Polygon(n, vertice)
    
# Optimized with Memoization
class Triangulation:
    def __init__(self):
        self.table = {}

    def triangulate(self, poly):
        """ Find all possible triangulations """
            # This is already a triangle
        if poly in self.table:
            return self.table[poly]
        
        if poly.n == 3:
            self.table[poly] = [[]]
            return [[]]
        
        output = []
        # This the edge from index 0 to index 1
        fixed_base = poly.edges[0]
        
        for i in range(2, poly.n):
            current_vertex = poly.vertice[i]
            line_1 = (fixed_base[0], current_vertex)
            line_2 = (fixed_base[1], current_vertex)
            
            if is_line_segment_in_polygon(line_1, poly) and is_line_segment_in_polygon(line_2, poly):
                # This means that this is a valid triangle!
                if i == 2:
                    diagonals = [line_1]
                elif i == poly.n - 1:
                    diagonals = [line_2]
                else:
                    diagonals = [line_1, line_2]
                
                # Create the two polygons divided by this triangle
                left_index = list(range(1, i + 1))
                right_index = list(range(i, poly.n)) + [0]
                
                left_triangulations = []
                right_triangulations = []
                
                # This means a left polygon exists
                if len(left_index) > 2:
                    left_vertices = list(map(lambda i: poly.vertice[i], left_index))
                    left_polygon = Polygon(len(left_index), left_vertices)
                    left_triangulations = self.triangulate(left_polygon)
                
                # This means a right polygon exists
                if len(right_index) > 2:
                    right_vertices = list(map(lambda i: poly.vertice[i], right_index))
                    right_polygon = Polygon(len(right_index), right_vertices)
                    right_triangulations = self.triangulate(right_polygon)
                
                if left_triangulations != [] and right_triangulations != []:
                    # Combining the Triangulations
                    for left_diagonals in left_triangulations:
                        for right_diagonals in right_triangulations:
                            possible_triangulation = left_diagonals + diagonals + right_diagonals
                            output.append(possible_triangulation)
                elif left_triangulations != [] and right_triangulations == []:
                    for left_diagonals in left_triangulations:
                        possible_triangulation = left_diagonals + diagonals
                        output.append(possible_triangulation)
                elif right_triangulations != [] and left_triangulations == []:
                    for right_diagonals in right_triangulations:
                        possible_triangulation = diagonals + right_diagonals
                        output.append(possible_triangulation)
                else:
                    # Shouldn't reach this case
                    print("Shouldn't reach this case")
        
        self.table[poly] = output               
        return output

if __name__ == "__main__":
    # square = Polygon(4, [Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)])
    # print(is_line_segment_in_polygon([Point(0, 0,), Point(1, 1.1)], square))
    # polygon = Polygon(5, [Point(0, 0), Point(2, 0), Point(2, 2), Point(1, 1), Point(0, 2)])
    polygon = regular_polygon(6)
    visualize_polygon(polygon)
    trig = Triangulation()
    
    iter = 1
    for i in range(iter):
        print("Iteration ", i)
        start = time.time()
        output = trig.triangulate(polygon)
        # output = find_all_triangulations(polygon)
        print("Number of Triangulations: ", len(output))
        end = time.time()
        print("Time: ", end - start)
        
    for i, o in enumerate(output):
        visualize_triangulation(polygon, o, "Polygon Number " + str(i))