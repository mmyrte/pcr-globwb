# README: Documentation to the adjustments to the snow and glacier modules.
### Joren Janzing (joren.janzing@slf.ch)
### Version 
18.11.2025
### Related to:
Janzing, J., Wanders, N., van Tiel, M., van Jaarsveld, B., Karger, D. N., & Brunner, M. I. (2024). Hyper-resolution large-scale hydrological modelling benefits from improved process representation in mountain regions. EGUsphere, 2024, 1-46.

## General
We have made several adjustments to the snow module used in PCR-GLOBWB and added a new glacier module as well.\
We have added two new files:
- `glaciers.py` 
- `snow.py`

We also have made some additional changes so that these files are called from the landSurface.py and landCover.py files.

For more detailed information, see: \
Janzing, J., Wanders, N., van Tiel, M., van Jaarsveld, B., Karger, D. N., & Brunner, M. I. (2024). \
Hyper-resolution large-scale hydrological modelling benefits from improved process representation in mountain regions. EGUsphere, 2024, 1-46.

## Prerequisites
All of this runs within the context of the 30-arcsec PCR-GLOBWB model, which is a global hydrological model developed at Utrecht University.\
It has the same dependencies as the original model, including among others the pcraster library.

## Contents:
### `snow.py`
#### `initializeSnow`
    Initialize the snow module: reads boolean variables from configuration file and translates these to the code. If snow redistribution is applied, it can also load additional relevant data.

#### `updateSnowFall`
    Divide the incoming precipitation into snowfall and rainfall. This uses the snow-to-rain transition range from Magnusson et al., 2014.

#### `updateSnowFall`
    Divide the incoming precipitation into snowfall and rainfall. This uses the snow-to-rain transition range from Magnusson et al., 2014.
   
#### `snowMeltSlaterAndClark`
    Running the actual snow module. Based on the initial condition specified in the configuration files, it can make the degree-day factor (DDF) vary with season (Slater and Clark, 2006), or albedo and add a precipitation melt term (Kraaijenbrink et al., 2021).

#### `simplifiedFreyAndHolzmann_pcraster`
    Runs a simple snow redistribution scheme based on Frey and Holzmann, 2015.

### `glaciers.py`
#### `initializeGlacier`
    Initialize the glacier module: reads boolean variables from configuration file and translates these to the code. Depending on the glacier module chosen, it also reads additional relevant input data.

#### `updateStaticGlacier`
    Day to day updates of the glacier module: calculating glacier melt and accumulation.

#### `glacierSlideImmerzeel`
    Laterally transporting the glacier ice following Immerzeel et al., 2012.

#### `glacierSlideImmerzeel`
    Laterally transporting the glacier ice following Immerzeel et al., 2012.

#### `updateDeltaH`
    Updating the delta H approach at the end of the year based on Huss et al, 2010 and Seibert et al., 2018.

#### `initializeDeltaH`
    Initializing the Delta-H approach and defining the glaciers shapes when a specific fraction of the glaciers has lost its mass based on Huss et al, 2010 and Seibert et al., 2018.


## References
Frey, S., & Holzmann, H. (2015). A conceptual, distributed snow redistribution model. Hydrology and Earth System Sciences, 19(11), 4517-4530.\
Huss, M., Jouvet, G., Farinotti, D., & Bauder, A. (2010). Future high-mountain hydrology: a new parameterization of glacier retreat. Hydrology and Earth System Sciences, 14(5), 815-829.\
Immerzeel, W. W., Van Beek, L. P. H., Konz, M., Shrestha, A. B., & Bierkens, M. F. P. (2012). Hydrological response to climate change in a glacierized catchment in the Himalayas. Climatic change, 110(3), 721-736.\
Janzing, J., Wanders, N., van Tiel, M., van Jaarsveld, B., Karger, D. N., & Brunner, M. I. (2024). Hyper-resolution large-scale hydrological modelling benefits from improved process representation in mountain regions. EGUsphere, 2024, 1-46.\
Kraaijenbrink, P. D., Stigter, E. E., Yao, T., & Immerzeel, W. W. (2021). Climate change decisive for Asia’s snow meltwater supply. Nature Climate Change, 11(7), 591-597.\
Magnusson, J., Gustafsson, D., Hüsler, F., & Jonas, T. (2014). Assimilation of point SWE data into a distributed snow cover model comparing two contrasting methods. Water resources research, 50(10), 7816-7835.\
Seibert, J., Vis, M. J., Kohn, I., Weiler, M., & Stahl, K. (2018). Representing glacier geometry changes in a semi-distributed hydrological model. Hydrology and Earth System Sciences, 22(4), 2211-2224.\
Slater, A. G., & Clark, M. P. (2006). Snow data assimilation via an ensemble Kalman filter. Journal of Hydrometeorology, 7(3), 478-493.\