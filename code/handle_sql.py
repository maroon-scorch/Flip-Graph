import sqlite3
import pandas as pd
from polygon import Polygon, Point, visualize_polygon, Triangulation, visualize_triangulation
from flip_graph import flip_graph, visualize_graph

sql_link = 'graph.db'

class SQLHandler:
    def __init__(self, link):
        self.link = link
        # Polygon is represented as string literal of its points
        # Graph is the adjancey matrix written in JSON literal
        self.header = [
            ["polygon", "TEXT"],
            ["vertice", "INT"],
            ["is_convex", "NUMBER(1)"],
            ["flip_graph", "TEXT"]
        ]
        self.keys = [tuple[0] for tuple in self.header]
        
        self.conn = sqlite3.connect(self.link)
    
    def initialize_table(self):      
        c = self.conn.cursor()  
        header_row = [tuple[0] + " " + tuple[1] for tuple in self.header]
        row = ''
        for i in range(len(header_row)):
            if i != len(header_row) - 1:
                row += header_row[i] + ", "
            else:
                row += header_row[i]
        print(row)
        
        table_command = 'CREATE TABLE IF NOT EXISTS graph (' + row + ')'
        c.execute(table_command)
        self.conn.commit()
        
    def add_entry(self, polygon, graph):
        c = self.conn.cursor()
                
        row = pd.DataFrame({
            'polygon': [str(polygon)],
            'vertice': [len(polygon.vertice)],
            'is_convex': [polygon.is_convex()],
            'flip_graph': [str(graph)]
        })
        
        row.to_sql("graph", self.conn,if_exists='append',index=False)
        # c.execute("INSERT INTO graph (polygon,vertice,is_convex,flip_graph) values(?,?,?,?)", row['polygon'], row['vertice'], row['is_convex'], row['flip_graph'])
        self.conn.commit()

if __name__ == "__main__":
    print("--- Initalizing the Data Base ---")
    test = SQLHandler(sql_link)
    # test.initialize_table()
       
    polygon = Polygon(5, [Point(0, 0), Point(2, 0), Point(2, 2), Point(1, 1), Point(0, 2)])
    visualize_polygon(polygon)
    graph = flip_graph(polygon)
    test.add_entry(polygon, graph)
    
 