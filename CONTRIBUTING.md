# Contributing to APackOfTheClones Development

The goal of this guide is to help you get up and contributing to APackOfTheClones as quickly as possible. Note that this guide is a modified version of ggplot2's github `contributing.md`. The guide is divided into two main pieces:

1. Filing a bug report in an issue.
2. Suggesting a change via a pull request.

Please note that ggplot2 is released with a [Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this project, 
you agree to abide by its terms.

## Bug Reports

When filing a bug report via [github issues](https://github.com/Qile0317/APackOfTheClones/issues), the most important thing is to include a minimal 
reproducible example if possible, so that the problem can be verified quickly, and then figure out how to fix it. There are three things you need to include to make your example reproducible: required packages, data, code.

1.  **Packages** should be loaded at the top of the script, so it's easy to
    see which ones the example needs.
  
2.  The easiest way to include **data** is to use `dput()` to generate the R code 
    to recreate it. For example, to recreate the `mtcars` dataset in R,
    I'd perform the following steps:
  
       1. Run `dput(mtcars)` in R
       2. Copy the output
       3. In my reproducible script, type `mtcars <- ` then paste.
       
    However, as this package intends to work with Seurat objects that tend to get massive, and the data privacy may be an issue, it is recommended to try and reproduce the bug on a potentially modified version of the example seurat object included in this package with ```data("combined_pbmc")```.
  
3.  Spend a little bit of time ensuring that your **code** is easy for others to
    read:
  
    - make sure you've used spaces and your variable names are concise, but
      informative
  
    - use comments to indicate where your problem lies
  
    - do your best to remove everything that is not related to the problem.  
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a 
fresh R session and pasting your script in.

(Unless you've been specifically asked for it, please don't include the output 
of `sessionInfo()`. Although if you feel the issue is for sure related to versioning, go ahead.)

## Pure Feature Requests

Use github issues to file pure enhancements/feature requests, and describe in as much detail as possible what functionality/optimizations could be made for which functions, and potentially pseudocode or actual code for the implementation. Additionally motivate with why it would be a positive contribution for this package.

## Pull requests

To contribute a change, you follow these steps:

1. Create a branch in git and make your changes.
2. Push branch to github and issue pull request (PR).
3. Discuss the pull request.
4. Iterate until either the PR is accepted or decided that it's not
   a good fit.

Each of these steps are described in more detail below. This might feel 
overwhelming the first time you get set up, but it gets easier with practice. 
If you get stuck at any point, please reach out for help to the contributors.

If you're not familiar with git or github, please start by reading <http://r-pkgs.had.co.nz/git.html>

Pull requests will be evaluated against a seven point checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivate the
    need for change.

    Also include this motivation in `NEWS` so that when a new release of
    ggplot2 comes out it's easy for users to see what's changed. Add your
    item at the top of the file and use markdown for formatting. The
    news item should end with `(@yourGithubUsername, #the_issue_number)`.

2.  __Only related changes__. Before you submit your pull request, please
    check to make sure that you haven't accidentally included any unrelated
    changes. These make it harder to see exactly what's changed, and to
    evaluate any unexpected side effects.

    Each PR corresponds to a git branch, so if you expect to submit
    multiple changes make sure to create multiple branches. If you have
    multiple changes that depend on each other, start with the first one
    and don't submit any others until the first one has been processed.

3.  __Use the coding style of this package__. Maintaining a consistent style across the whole code base makes it much easier to jump into the code. If you're modifying existing code, a separate pull request to change the style would be greatly appreciated. In general, this would be a breaking change if it is an externally exported function, which ideally should be avoided.

4.  If you're adding new parameters or a new function, you'll also need
    to document them with [roxygen](https://github.com/klutometis/roxygen).
    Make sure to re-run `devtools::document()` on the code before submitting.

5.  If fixing a bug or adding a new feature to a non-graphical function,
    please add a [testthat](https://github.com/r-lib/testthat) unit test.

6.  If fixing a bug in the visual output, please add a visual test with
    ```expect_doppelganger()```, but try as much as possible to reduce
    the svg snapshot size by minimizing the number of individual polygons
    displayed.

7.  If you're adding a new graphical feature, please add a short example
    to the appropriate function.

This seems like a lot of work but don't worry if your pull request isn't perfect.
A pull request ("PR") is a process, that require review and approval before merging.
