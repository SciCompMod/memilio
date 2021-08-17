## Merge Request - GuideLine Checklist 

**Guideline** to check code before resolve WIP and approval, respectively.
As many checkboxes as possible should be ticked.

### Checks by code author:
* [ ] There is at least one issue associated with the pull request.
* [ ] The branch follows the naming conventions as defined in the [git workflow](git-workflow).
* [ ] New code adheres with the [coding guidelines](coding-guidelines)
* [ ] Tests for new functionality has been added
* [ ] A local test was succesful
* [ ] There is appropriate **documentation** of your work. (use doxygen style comments)
* [ ] If new third party software is used, did you pay attention to its license? Please remember to add it to the wiki after successful merging.
* [ ] If new mathematical methods or epidemiological terms are used, has the glossary been updated ? Did you provide further documentation ?
 is present or referenced. Please provide your references.
* [ ] The following questions are addressed in the documentation*:  Developers (what did you do?, how can it be maintained?), For users (how to use your work?), For admins (how to install and configure your work?)
* For documentation: Please write or update the Readme in the current working directory!

### Checks by code reviewer(s):
* [ ] Is the code clean of development artifacts e.g., unnecessary comments, prints, ...
* [ ] The ticket goals for each associated issue are reached or problems are clearly addressed (i.e., a new issue was introduced).
* [ ] There are appropriate **unit tests** and they pass. The meaning of "appropriate" is defined by the [test strategy](not yet defined).
* [ ] The git history is clean and linearized for the merge request.
* [ ] Coverage report for new code is acceptable. 

