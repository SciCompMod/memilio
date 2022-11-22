## Merge Request - GuideLine Checklist 

**Guideline** to check code before resolve WIP and approval, respectively.
As many checkboxes as possible should be ticked.

### Checks by code author:
Always to be checked:
* [ ] There is at least one issue associated with the pull request.
* [ ] The branch follows the naming conventions as defined in the [git workflow](git-workflow).
* [ ] New code adheres with the [coding guidelines](coding-guidelines)
* [ ] Make sure that the [pre-commit linting/style checks pass](https://github.com/DLR-SC/memilio/wiki).

If functions were changed or functionality was added:
* [ ] Tests for new functionality has been added
* [ ] A local test was succesful

If new functionality was added:
* [ ] There is appropriate **documentation** of your work. (use doxygen style comments)

If new third party software is used:
* [ ] Did you pay attention to its license? Please remember to add it to the wiki after successful merging.

If new mathematical methods or epidemiological terms are used:
* [ ] Are new methods referenced? Did you provide further documentation? Has the glossary been updated? 

The following questions are addressed in the documentation if need be: 
* [ ] Developers (what did you do?, how can it be maintained?)
* [ ] For users (how to use your work?)
* [ ] For admins (how to install and configure your work?)

* For documentation: Please write or update the Readme in the current working directory!

### Checks by code reviewer(s):
* [ ] Is the code clean of development artifacts e.g., unnecessary comments, prints, ...
* [ ] The ticket goals for each associated issue are reached or problems are clearly addressed (i.e., a new issue was introduced).
* [ ] There are appropriate **unit tests** and they pass.
* [ ] The git history is clean and linearized for the merge request. All reviewers should squash commits and write a simple and meaningful commit message.
* [ ] Coverage report for new code is acceptable. 

