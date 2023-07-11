Check [git workflow](https://github.com/DLR-SC/memilio/wiki/git-workflow) and use the **draft** feature before a Pull request/Merge request.

# Changes and Information

Please briefly list the changes made and what the Reviewer should look out for, as well as additional Information:

[Describe here in 2–3 lines. For details, reference to issues to avoid redundant information.]

## Checks by code author

- [ ] Every addressed issue is linked (use the "Closes #ISSUE" keyword below)
- [ ] New code adheres to [coding guidelines](https://github.com/DLR-SC/memilio/wiki/Coding-guidelines)
- [ ] No large data files have been added (files should in sum not exceed 100 KB, avoid PDFs, Word docs, etc.)
- [ ] Tests are added for new functionality and a local test run was successful
- [ ] Appropriate **documentation** for new functionality has been added (Doxygen in the code and Markdown files if necessary)
- [ ] Proper attention to licenses, especially no new third-party software with conflicting license has been added

## Checks by code reviewer(s)

- [ ] Corresponding issue(s) is/are linked and addressed
- [ ] Code is clean of development artifacts (no deactivated or commented code lines, no debugging printouts, etc.)
- [ ] Appropriate **unit tests** have been added, CI passes and code coverage is acceptable (did not decrease)
- [ ] No large data files added in the whole history of commits(files should in sum not exceed 100 KB, avoid PDFs, Word docs, etc.)
