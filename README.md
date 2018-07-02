# usrail
A demonstration of using minimum-spanning trees and traveling-salesman solutions to design a national rail network for the United States 

## Requirements

Python 2.7

### Required Python Packages
  * [NumPy](http://www.numpy.org/)
  * [SciPy](https://www.scipy.org/)
  * [Matplotlib](https://matplotlib.org/)
  * [pyshp](https://github.com/GeospatialPython/pyshp)
  * [Pandas](https://pandas.pydata.org/)
  * [PyConcorde](https://github.com/jvkersch/pyconcorde)
 
 ## Motivation
 
 This was inspired by some posts on the New Urbanist Memes for Transit-Oriented Teens Facebook group about ideal Amtrak configurations. Having had some practice with optimization and network tools, I decided to see what would happen if I let these algorithms loose on the geography of the United States.
 
 ## Problem Definition
 
 Optimizing a network of vehicles in reality is prohibitively difficult. Exploring the problem requires some simiplifying assumptions. To simplify the problem, I attempt to answer the question: 
 
 What is the minimum configuration of a rail network to connect all 378 (metropolitan statistical areas](https://en.wikipedia.org/wiki/Metropolitan_statistical_area) (MSAs) in the contiguous United States that
 * minimizes the total amount of track required,
 * maximizes the number of passengers on direct intercity connections, *or* 
 * maximizes revenue generated from passengers on direct intercity connections?
 
 The scope of the problem is defined so it can be solved by finding the  [minimum spanning tree](https://en.wikipedia.org/wiki/Minimum_spanning_tree) (MST) of the MSAs with edge weights defined by distance, passenger volume, and revenue passenger-distance. The key properties of the MST solution are:
 
  * The sum of the edge weights is minimized or maximized  
 * There is exactly one path between any two cities
 * There are no loops
 * The sum of the edge weights is minimized or maximized
 
 It's important to keep in mind that the MSTs represent the minimum sufficient solution to this problem. Introducing a connection between two close cities that aren't on the same branch of the tree is possible, but outside of the scope of this problem as defined.
 
 ## Method and Results
 
 ### Data
 
The network nodes are based on the geometric centroids of all 378 metropolitan statistical areas of the continuguous United States. This excludes 21 MSAs in Alaska, Hawaii, and Puerto Rico. The MSA boundaries come from the Census Bureau and are included in this repo. The MSAs I used contained 86% of the US population in 2017. 

 ### Weighting Schemes
 
For all solutions, I use a simple gravity model to represent intercity traffic volume. The model is based on the understanding that the activity between two cities is proportional to the products of their respective populations and inversely proportional to the distance between them. The larger any pair of cities is, the more people move between them, but the amount of traffic between two cities should fall with the distance between them. 
 
 The gravity model I use is just 
 
 ```
        p_i*p_j
 w_ij ~ ———————,
         d_ij^2
 ```
 where ```w``` is the intercity traffic volume between cities ```i``` and ```j```, or the "edge weight" in graph theory parlance, ```d_ij``` is the distance between the two cities, and ```p_i``` and ```p_j``` are the cities' respective populations.
 
 The gravity model is a limiting simplification. I also don't adjust the intercity traffic based on the final configuration of the network. In reality, the creation of a rail connection between cities could increase the volume of travel between them. I also don't consider the effects the length of a route or number of connections might have on the attractiveness of the route to passengers.
 
 #### Minimum Track Length
 
 Here, the metric is the distance between the two cities. It may be understood as a reflection of the cost of constructing and maintaining a rail network. Minimizing the total track used in connecting every city in the United States is the most intuitive solution I can generate.
 
  #### Minimum Track Length Loop
  
  This is a closed loop that connects all cities on while using the minimum possible amount of track. This is a solution to the traveling salesman problem. I like to think of it as the "Snowpiercer" solution. I throw this one in for fun.
 
 #### Maximum Ridership
 In this solution, the metric is the number of passengers on direct connections between cities. This metric may also be thought of as the revenue per unit track length, a reflection of the profitability of a connection. 
 
  #### Maximum Revenue
  Revenue in transportation is proportional to the amount of material moved between two points and the distance between those points. In this weighting scheme, I simply multiply the gravity model by the distance between cities. Because construction and operating costs and are proportional to total track mileage and I do not have a microeconomic model built into this solution, this metric may not reflect the actual profitability of a connection. 
  
  ### Results
 
 In the solutions below, the populations of each city are shown in yellow circles sized according to population. Connections between cities are represented by solid black lines. MSA boundaries are shown in dark gray. Micropolitan areas are represented in light gray for reference, but do not factor into the solutions.
 
 ![Solutions](https://github.com/ryanahardy/usrail/blob/master/solutions.png "Rail Network Solutions")
 
  When reporting total ridership and revenue in the optimization results, I refer to normalized scores instead of dimensionally meaningful units. These are obtained by dividing each result by the mean of all results.

 ![Scores](https://github.com/ryanahardy/usrail/blob/master/scores.png "Rail Network Scores")
 
 The minimum track length solution yields a network with 35,700 km of track. It has a ridership score of 0.7 and a revenue score of 0.4. Amtrak, by comparison, uses 34,400 km of track.
 
 The minimum track length loop (traveling salesman or Snowpiercer solution) requires a slightly larger 41,500 km of track. Its respective ridership and revenue scores of 0.7 and 0.4 are comparable to its MST cousin.
 
 Maximizing ridership prioritizes connections between major cities. In this solution, a hub-and-spoke model emerges. This solution requires 79,000 km of track and has a ridership score of 1.4 and a revenue score of 1.5, twice what is offered by minimizing track length.
 
 Maximizing revenue rewards connections between large cities across large distances. Consequently, the total track required is 162,000 km, over four times the minimum amount required. The ridership score drops to 1.2 relative to the maximum ridership solution, but the revenue score increases to 1.7. At this point, the solution is a pure hub-and-spoke model. The solution generates connections that cross open ocean, making them ideally suited for aircraft rather than trains. Unlike other solutions shown here, most of the direct connections are used by exactly one line, reducing the profitability of its span of track. This is a problem airlines don't have to deal with, making the hub-and-spoke model ideal for air travel.
 
 An interesting feature of both hub-and-spoke models is how they divide the United States into distinct clusters. In the maximum revenue solution, the main New York hub connects to six secondary hubs in Chicago, Atlanta, Dallas—Forth Worth, Los Angeles, Detroit, and Miami. These hubs connect to tertiary hubs in places like Seattle, Houston, Denver, and San Francisco, which serve even smaller cities. One could draw boundaries around the secondary hubs and their spokes to define geographically distinct regions of the United States. That's a different project, however.
 

