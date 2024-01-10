
## R CMD check results

0 errors | 0 warnings | 1 note

## Additional notes

- The test failures from before have been diagnosed asfloating point precision issues leading to numeric comparisons producing different function outputs on CRAN, whereas they always worked on my local machine and CI platforms. The issue has now been fixed.
