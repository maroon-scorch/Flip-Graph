from tkinter import *
from tkinter import ttk
import time
from polygon import Point, Polygon, visualize_polygon
from flip_graph import flip_graph, visualize_graph

win=Tk()

point_list = []

screen_width = 400
screen_height = 400
ratio = 5

# Setup of the window
win.geometry("400x500")
win.resizable(False, False)
win.title('GUI')

prev_point = [None, None]

def draw_line(event):
   x1=event.x
   y1=event.y
   x2=event.x
   y2=event.y
   
   # Draw an oval in the given co-ordinates
   canvas.create_oval(x1,y1,x2,y2,fill="black", width=10)
   
   if prev_point != [None, None]:
       canvas.create_line(prev_point[0], prev_point[1], x1, y1, fill="black", width=5)

   prev_point[0] = x1
   prev_point[1] = y1
   cartX = x1 - screen_width/2
   cartY = screen_height/2 - y1
   cartPoint = Point(cartX/ratio, cartY/ratio)
   
   print(cartPoint)
   
   point_list.append(cartPoint)

canvas=Canvas(win, width=screen_width, height=screen_height, background="white")
canvas.grid(row=0, column=0)
canvas.bind('<Button-1>', draw_line)
click_num=0

def clickClearButton():
    point_list.clear()
    prev_point[0] = None
    prev_point[1] = None
    canvas.delete('all')
    
def clickRunButton():
    print(point_list)
    
    file = open("input.txt",'wt')
    file.write(str(len(point_list)) + "\n")
    
    for point in point_list:
        px = point.x
        py = point.y
        line = str(px) + " " + str(py) + "\n"
        file.write(line)
    file.close()
    
    poly = Polygon(len(point_list), point_list)
    win.destroy()
    visualize_polygon(poly)
    
    start = time.time()
    graph = flip_graph(poly)
    end = time.time()
    print("Time: ", end - start)
    
    
    file = open("output.txt",'wt')
    for i in range(len(graph.keys())):
        line = str(i) + " : " + str(graph[i])  + "\n"
        file.write(line)
    file.close()
    visualize_graph(graph)
    

clearButton = Button(text="Clear", command= clickClearButton)
clearButton.place(x=0, y=450)

runButton = Button(text="Run", command= clickRunButton)
runButton.place(x=200, y=450)

win.mainloop()