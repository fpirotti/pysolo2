The problem to solve is to select the optimal cities were to collect residual biomass from specific types of crops: vineyards and olive orchards.

First of all we apply two conditional filters with a binary threshold. A city is

$$
y =\begin{cases}1 & \text{if } x_1 \geq T_1 \;\text{and}\; x_2 \geq T_2 \\0 & \text{otherwise}\end{cases}
$$

Choosing a specific city implies specific advantages and disadvantages, (also referred to as gains and losses) when referring to an optimization strategy. The following factors have an impact which has a spatial connotation:

-   **D~cb~** - distance of the city to the nearest biomass plant: having a biomass plant close to the city is actually a negative aspect as biomass plants compete for biomass resources. - thus a direct positive correlation of suitability of the city with distance from existing biomass plants.
-   **D~cr~** - distance of the city to the nearest refinery: refineries imply an infrastructure for distribution of energy resources, thus the closer an existing refinery is, the better it is for the city. Therefore there is an inverse correlation with distance, i.e. less distance with a refinery is better.
-   **D~cc~** -distance of the city to the nearest company: the same as above, having a company close-by is a positive asset, therefore we have an inverse correlation with distance - i.e. the more distant an existing compay is the less suitable that city is.
-   **D~cp~** -distance of the city to the nearest ports: the same as the two point just above. A nearby port is a positive asset of the city as transportation to other users of bioenergy is easier.

All the above provide an initial "weight" to each cities.

-   

We first proceed to defining the location of vineyards and olive orchards using CORINE Land Cover at 100 m from 2018. The new High Resolution Land Cover Plus at 10 m resolution could potentially be used, but it should be noted that at 100 m Spain will have about 3.3E6 cells with such land cover, therefore at 10 m we would expect 100x more, making the computation more complex. This can be addressed in the future, for now we will test the method with 100 m pixels.

We can pose the problem considering that two things have to be defined: how many cities to pick from the total and which cities to pick. We will define the optimal decision as a function of gain and loss. *Gain* is how much residual biomass **g** can be collected in a city from all neighbouring cells with vineyards and olive orchards, and the *loss* consists in the economic cost of building the pysolo plant, which we can call "unit cost" **UC~pp~** (economic), the distance of each residual biomass cell from the nearest road infrastructure (**D~gr~**) to consider costs, the distance of the city to the nearest port (**D~cp~**) , the distance of the city to the nearest biorefinery (**D~cb~**) and the distance of the city to the nearest company (**D~cc~**). If we consider **k** the total number of cities we have, then we want to build the plant in **n** cities (**n \<= k**). The number of cities were to build the pysolo plant, **n,** depends on the loss/gain ratio from each possible combination. The loss factors are all constant except **UC~pp~** which we can modify by increasing it to see which cities are optimal in different **UC~pp~** scenarios. We potentially can test the following number of combinations **C** .

$$C(n,k) = \binom{n}{k} = \frac{n!}{k!(n-k)!}$$

Because to test all possible values of $n$ we would start from $k$ (all cities) to only choosing a single city, \$ n={k...1} \$, that would be computationally expensive.

To lower the computational effort we take another approach.

We start by assigning zero cost to building a pysolo plant, **UC~pp~** = 0, for the sake of starting the computation by "reductio ad absurdum". Then all cities will be used, except the ones that are not near any cell with residual biomass from vineyards or olive groves as they will have zero cost but also zero gain, so no marginal value. Of course zero cost is not realistic, so we start increasing the value of **UC~pp~** and check which cities to remove to have the best combination for our increasting costs.

The optimization function maximizes the combination of weighted benefits and costs by finding k cities by solving the following function:

${maximize}\quad \sum _{i=1}^{n}(g_{i}-l_{i})\cdot x_{i}$

$x_{i}\in \{0,1\}\quad \forall i\in \{1,\dots ,n\}$

The binary nature of the decision variables x is simply that the $n$ cities will get zero (0) or one (1) depending on if they are included as candidate cities or not, thus $k$ will be the number of cities with $1$. The system tries all combinations and finds the optimal one.

$g_{i}$: The *gain* (or value) associated with city $i$.

$$
g_i = f(sum(pixel suitability * 1/distanza) )
$$

$l_{i}$: The loss associated with selecting city $i$ (it is composed of cost of building plus other factors noted above)

$x_{i}$: A binary variable for each city equal to 1 if city $i$ is selected, and 0 otherwise.

$i$. ith city

$n$: The total number of cities to choose from.

Note 1 - the gain and the part of the loss "*Dgr*" for each city is calculated from the sum of the normalized NDVI values inversely weighted with the cell distance from the road infrastructure.

Note 2 - the normalized values scale was in the range of 0-255 in integers to improve memory efficiency without significantly loosing precision.

Note 3 - the limitation is that the gain and the part of loss determined by distance from roads of cells is assigned to each cell, and each cell is assigned to the nearest city using the Voronoy spatialization. It is reasonable to think that in the case of two or three neighbouring cities close to each other where only one is selected, the cells assigned to the cities that were not selected would be re-assigned to the select city. This is not considered in this procedure, due to the paradox that we know if cities are selected or not only after the optimization function is solved. We would have to re-assign the values of the gain due to more cells being added to the selected cities. We can reply that adding cells after the optimization would only potentially increase the gain, thus reinforce the decision.
