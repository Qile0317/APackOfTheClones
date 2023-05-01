# APackOfTheClones 0.1.0 (2023-05-01)
Initial semi-stable release. The main functions of the package are working with their default parameters. 

Here is a list of currently known bugs:
- Automated cluster repulsion uses a highly flawed mathematical formula loosely based on force directed graph drawing techniques, the resulting plots are not presentable and the function has been disabled temporarily
- The functions to add a legend of circle sizes do not work as intented and is temporarily disabled.\

Upcoming features are:
- Customizable cluster coloration, not just based on the original Seurat/ggplot palette
- User-controlled cluster shifting
- Comprehensive cluster repulsion
- legend placement
- automated optimization of the initial parameters of `clonal_expansion_plot`, especially the `clone_scale_factor`
- BCR library integration