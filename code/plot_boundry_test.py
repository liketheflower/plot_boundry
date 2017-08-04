from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import numpy as np


points = np.random.rand(300, 2)   # 30 random points in 2-D
hull = ConvexHull(points)

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(points[:,0], points[:,1], 'o')
for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
#plt.show()



delta = 0.4
points_list = points.tolist()
#rint points
length = len(points_list)
for i in range(length):
    points_list.extend([[points_list[i][0]-delta,points_list[i][1]-delta], [points_list[i][0]-delta,points_list[i][1]+delta],[points_list[i][0]+delta,points_list[i][1]-delta],[points_list[i][0]+delta,points_list[i][1]+delta]])
points =  np.asarray(points_list)
#plt.figure(2)
plt.plot(points[:,0], points[:,1], 'o')
hull = ConvexHull(points)
for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'r-')

plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
plt.show()


def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0
    
    
def plot_in_hull(p, hull):
    """
    plot relative to `in_hull` for 2d data
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection

    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    # plot triangulation
    poly = PolyCollection(hull.points[hull.vertices], facecolors='w', edgecolors='b')
    plt.clf()
    plt.title('in hull')
    plt.gca().add_collection(poly)
    plt.plot(hull.points[:,0], hull.points[:,1], 'o', hold=1)


    # plot the convex hull
    edges = set()
    edge_points = []

    def add_edge(i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(hull.points[ [i, j] ])

    for ia, ib in hull.convex_hull:
        add_edge(ia, ib)

    lines = LineCollection(edge_points, color='g')
    plt.gca().add_collection(lines)
       

    # plot tested points `p` - black are inside hull, red outside
    inside = in_hull(p,hull)
    print type(inside)
    print inside
    print np.logical_not(inside)
    
    plt.plot(p[ inside,0],p[ inside,1],'.k')
    plt.plot(p[np.logical_not(inside),0],p[np.logical_not(inside),1],'.r')  
    plt.show() 
    
    
tested = np.random.rand(1000,2)
cloud  = np.random.rand(50,2)

print in_hull(tested,cloud)

plot_in_hull(tested,cloud)

