# SteadyCom
(The functions in this repository are already incorporated into the **[COBRA toolbox](https://github.com/opencobra/cobratoolbox)** since 2018. It is recommended to install COBRA toolbox to use them. Check out this [tutorial](https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialSteadyCom.html).)

Matlab implementation for computing the maximum growth rate of a microbial community at community steady-state.

Chan SHJ, Simons MN, Maranas CD (2017) SteadyCom: Predicting microbial abundances while ensuring community stability. PLoS Comput Biol 13(5): e1005539. https://doi.org/10.1371/journal.pcbi.1005539

Main functions:

`SteadyCom`: given a community model constructed by createCommModel, find the maximum community growth rate by the SteadyCom procedure

`SteadyComFVA`: perform FVA under the SteadyCom framework

`SteadyComPOA`: perform pairwise FVA under the SteadyCom framework

Remark:
- Using the `ibm_cplex` solver for the COBRA toolbox is usally the fastest when running this code because the code directly calls the cplex object in Matlab
- Use `createMultipleSpeciesModel` to create a community model from individual models
- Use `getMultiSpeciesModelId` to retreive the IDs needed for `SteadyCom` from the created model
- See the [tutorial] for a detailed example
