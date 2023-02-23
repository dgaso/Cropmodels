# Soybean crop growth model

This code belongs to the paper: Predicting within-field soybean yield variability by coupling Sentinel-2 leaf area index with a crop growth model (https://doi.org/10.1016/j.agrformet.2021.108553)  

![image](https://user-images.githubusercontent.com/75320488/220874024-c8895541-c567-41db-a97d-a657061218a9.png)


The soybean model includes a crop and soil module. Soybean above-ground biomass is described with a daily time step over a unit field area of a square meter. The model uses a combination of water-driven productivity and radiation use efficiency to calculate the daily above-ground biomass increase. Dry matter accumulation in leaves, stems and seeds was defined using the static partitioning tables approach from WOFOST. Green leaf area index progress between emergence and beginning pod was obtained by multiplying leaf biomass produced each day and the specific leaf area (SLA) parameter for soybean. Yield formation starts slowly after the flowering stage until maturity. Two sources for pod growth are considered for computing biomass allocated in seed organs. The first source is directly derived by multiplying the growth rate after flowering time and the fraction partitioned to seed. 

The soil module is the simplified water balance described in: Campbell, G.S., Diaz, R., 1988. Simplified Soil-Water Balance Models to Predict Crop
Transpiration, in: Drought Research Priorities for the Dryland Tropics. pp. 15–26.
The soil profile is divided into N layers. Root exploration in the soil profile was described by a logistic curve and the fraction of roots in a layer depends on the root depth as rooting density decreases linearly with depth. Root water uptake (actual transpiration rate) is computed from each soil layer. Water uptake from each layer is assumed directly proportional to the difference in water potential between the soil in that layer and the xylem and is inversely proportional to the root resistance. Root resistance is proportional to the fraction of roots in a given layer. An estimate of the xylem water potential is computed by using water potential weighted by the fraction of roots in a given layer, potential transpiration rate and the total root resistance. The total root resistance decreases with plant growth, which is achieved by dividing the minimum resistance parameter by the proportion of intercepted light. The model does not simulate water movements and redistribution due to osmotic suction gradients. The first soil layer is only used to compute soil evaporation and to enter the rainfall into the soil profile, stored water in this layer is considered not to be available to the crop.


A detailed description of the model can be found in: https://doi.org/10.1016/j.agrformet.2021.108553 
