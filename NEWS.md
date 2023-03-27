# Pre-release notes
- tests unfinished
- github workflows unfinished

- In the future more features will be added, such as allowing BCR data integration as well

# APackOfTheClones 0.1.0 (2023-03-28)
Initial barebones release. The main functions of the package are working well with their default parameters. Here is a list of currently known bugs:
- Automated cluster repulsion uses a highly flawed mathematical formula loosely based on force directed graph drawing techniques. The user should ignore any repulsion related arguments. 
- The progress bars for individual cluster circle packing sometimes exceed 100% and/or are duplicated

Upcoming features are:
- Customizable cluster coloration, not just based on the original seurat/ggplot palette
- User-controlled cluster shifting
- Comprehensive cluster repulsion