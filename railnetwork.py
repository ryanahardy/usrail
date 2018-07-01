#!/usr/bin/env python
'''A demonstration of using minimum-spanning trees and traveling-salesman
solutions to design a national rail network for the United States.'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import shapefile
import pandas as pd
import scipy.sparse.csgraph as csg
from mpl_toolkits.basemap import Basemap
from concorde.tsp import TSPSolver

__author__ = "Ryan A. Hardy"
__copyright__ = "2018, Ryan A. Hardy"
__credits__ = ["Ryan A. Hardy"]
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Rob Knight"
__email__ = "hardy.r@gmail.com"
__status__ = "Development"

def plot_poly(m, sf, y=None, mask=None, cmap=plt.cm.RdBu_r, vmin=None, vmax=None, **kwargs):
    patches = []
    colors = []
    if y==None:
        y = np.ma.masked_invalid(np.zeros(len(sf.shapes()))+np.nan)
    if vmin==None:
        vmin = y.min()
    if vmax==None:
        vmax = y.max()
    for i in range(len(sf.shapes())):
        if mask is not None:
            if ~mask[i]:
                continue
        #print sf.records()[i][5]
        norm = (y[i]-vmin)/(vmax-vmin)
        shape = sf.shapes()[i]
        points = np.array(shape.points)
        if len(shape.parts) > 1:
            for j in range(len(shape.parts)):
                if j < len(shape.parts)-1:
                    poly = points[shape.parts[j]:shape.parts[j+1]]
                else:
                    poly = points[shape.parts[j]:]
                xy = np.array(m(poly[:, 0], poly[:, 1])).T
                patches.append(Polygon(xy, closed=False,
                                    facecolor=cmap(norm), **kwargs))
                plt.gca().add_patch(patches[-1])
                colors.append(y[i])
        else:
            poly = points
            xy = m(poly[:, 0], poly[:, 1])
            patches.append(Polygon(np.array(xy).T, closed=False,
                                facecolor=cmap(norm), **kwargs))
            plt.gca().add_patch(patches[-1])
            colors.append(y[i])
    p = PatchCollection(patches, cmap=cmap)
    p.set_array(colors)
    p.set_norm(plt.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=False))
    plt.gca().autoscale_view()
    return p



#sf_csa = shapefile.Reader("cb_2017_us_csa_5m/cb_2017_us_csa_5m")

#Load CBSA and US State shapefiles

sf_cbsa = shapefile.Reader("data/cb_2017_us_cbsa_20m/cb_2017_us_cbsa_20m")
sf_state = shapefile.Reader("data/cb_2017_us_state_20m/cb_2017_us_state_20m")
m = Basemap(projection='aea', lon_0=-96, lat_0=37.5, width=5e6, height=4e6)

#Isolate and compute MSA centroids

metro = (np.array(sf_cbsa.records())[:, -3]=='M1')
metro_where = np.where(metro)[0]
nmsa = metro.sum()

centroids = np.zeros((nmsa, 2))
x_c = np.zeros(nmsa)
y_c = np.zeros(nmsa)
for i in range(nmsa):
    points = np.array(sf_cbsa.shapes()[metro_where[i]].points)
    points_proj = np.array(m(*points.T)).T
    cross = np.diff(points_proj*np.roll(np.fliplr(points_proj), 1, axis=0), 1)
    A = cross.sum()*0.5
    x_c[i], y_c[i] = ((points_proj+np.roll(points_proj, 1, axis=0))*cross).sum(0)/(6*A)
centroids = np.array(m(x_c, y_c, inverse=True)).T

visible = (x_c < m.xmax)*(x_c > m.xmin)*(y_c < m.ymax)*(y_c > m.ymin)

centroids = centroids[visible]
nmsa = len(centroids)

#Load CBSA populations and match with shapefiles
csapop = pd.read_csv("data/PEP_2017_GCTPEPANNR.US23PR/PEP_2017_GCTPEPANNR.US23PR_with_ann.csv")
ref = np.where(np.matrix(np.array(sf_cbsa.records())[:, 3][metro][visible]).T==np.matrix(csapop['GC.target-geo-id2'].ix[2:].values))[1]
pop = np.array(csapop['respop72017'].ix[2:].values[ref], int)

#Compute intercity geographic distance matrix
dmat = np.zeros((nmsa, nmsa))
R = 6371. #km, radius of Earth
for i in range(nmsa):
    dmat[i] = R*np.arccos(np.sin(centroids[i, 1]*np.pi/180)*np.sin(centroids[:, 1]*np.pi/180)+np.cos(centroids[i, 1]*np.pi/180)*np.cos(centroids[:, 1]*np.pi/180)*np.cos((centroids[i, 0]-centroids[:, 0])*np.pi/180))
dmat[dmat!=dmat] = 0

#Define metrics based on gravity model
distance = dmat
ridership = np.nan_to_num(np.outer(pop, pop)/dmat**2) #Passengers
revenue = np.nan_to_num(ridership*distance) #Passenger-distance

#Compute minimum-spanning trees. Negate metrics for maximization.
names = ['Maximum Ridership', 'Maximum Revenue',
        'Minimum Track Length', 'Minimum Track Length\nLoop']
mst = [ csg.minimum_spanning_tree(-ridership),
        csg.minimum_spanning_tree(-revenue),
        csg.minimum_spanning_tree(distance)]

#Solve the traveling salesman problem
solver = TSPSolver.from_data(centroids[:, 1], centroids[:, 0], norm="GEO")
tour_data = solver.solve(verbose=False)

scores = pd.DataFrame(index=names, columns=['Ridership', 'Revenue', 'Track Length'])

#Report scores
for i in range(4):
    print names[i]
    if i < 3:
        I, J = np.where(mst[i].toarray()!=0)
    else:
        I =  [tour_data.tour[j] for j in range(nmsa)]
        J = np.roll(I, 1)
    print "Ridership: %f" % (ridership[I, J].sum()/1e10)
    print "Revenue: %f" % (revenue[I, J].sum()/1e12)
    print "Length: %f km" % distance[I, J].sum()
    print
    scores.iloc[i] = ridership[I, J].sum(), revenue[I, J].sum(), distance[I, J].sum()

#Normalize scores and plot
(scores/scores.mean()).iloc[[2, 3, 0, 1]].plot(kind='bar', rot=0, lw=0, figsize=(12, 6))
plt.ylabel("Normalized Score", weight='bold')
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
plt.grid('on', axis='y')
plt.legend(loc=0, frameon=False, ncol=3)
plt.xticks(np.arange(4), np.array(names)[[2, 3, 0, 1]], weight='bold')
plt.savefig("scores")


#Plot maps
plt.figure(figsize=(12.5, 10))
for j in range(3):
    plt.subplot(2, 2, [3, 4, 1][j])
    I, J = np.where(mst[j].toarray()!=0)
    plt.title(names[j], weight='bold', va='top')
    for i in range(nmsa-1):
        x, y = m.gcpoints(centroids[I[i], 0], centroids[I[i], 1], centroids[J[i], 0], centroids[J[i], 1], 10)
        m.plot(x, y, color='k')
    plot_poly(m, sf_state, fc='0.75', ec='1')
    plot_poly(m, sf_cbsa, mask=~metro, ec='0.75', fc='0.67', lw=0.5)
    plot_poly(m, sf_cbsa, mask=metro, ec='0.75', fc='0.5', lw=0.5)
    m.scatter(centroids[:, 0], centroids[:, 1], s= np.array(pop, int)/1e5, latlon=True, zorder=10, lw=0, c='lightyellow')

    plt.axis('off')

plt.subplot(2, 2, 2)
plt.title(names[3], weight='bold', va='top')
for i in range(nmsa-1):
    x, y = m.gcpoints(centroids[tour_data.tour[i], 0],
                        centroids[tour_data.tour[i], 1],
                        centroids[tour_data.tour[i+1], 0],
                        centroids[tour_data.tour[i+1], 1], 10)
    m.plot(x, y, color='k')
x, y = m.gcpoints(centroids[tour_data.tour[-1], 0],
                    centroids[tour_data.tour[-1], 1],
                    centroids[tour_data.tour[0], 0],
                    centroids[tour_data.tour[0], 1], 10)
m.plot(x, y, color='k')
plot_poly(m, sf_state, fc='0.75', ec='1')
plot_poly(m, sf_cbsa, mask=~metro, ec='0.75', fc='0.67', lw=0.5)
plot_poly(m, sf_cbsa, mask=metro, ec='0.75', fc='0.5', lw=0.5)
m.scatter(centroids[:, 0], centroids[:, 1], s= np.array(pop, int)/1e5, latlon=True, zorder=10, lw=0, c='lightyellow')
plt.axis('off')
plt.subplots_adjust(left=0, hspace=0, wspace=0, right=1, top=0.95, bottom=0)
plt.savefig("solutions")
plt.savefig("solutions", format='pdf')
