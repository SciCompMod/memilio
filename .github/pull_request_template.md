# Changes and Information

Please **briefly list the changes** (main added features, changed items, or corrected bugs) made:

-
-
-


If need be, add additional information and what the reviewer should look out for in particular:

-
-

## Merge Request - Guideline Checklist

Please check our [git workflow](https://memilio.readthedocs.io/en/latest/development.html#git-workflow). Use the **draft** feature if the Pull Request is not yet ready to review.

### Checks by code author

- [ ] Every addressed issue is linked (use the "Closes #ISSUE" keyword below)
- [ ] New code adheres to [coding guidelines](https://memilio.readthedocs.io/en/latest/development.html#coding-guidelines)
- [ ] No large data files have been added (files should in sum not exceed 100 KB, avoid PDFs, Word docs, etc.)
- [ ] Tests are added for new functionality and a local test run was successful (with and without OpenMP)
- [ ] Appropriate **documentation** for new functionality has been added (Doxygen in the code and explanations in the online documentation)
- [ ] Proper attention to licenses, especially no new third-party software with conflicting license has been added
- [ ] (For ABM development) Checked [benchmark results](https://memilio.readthedocs.io/en/latest/development.html#agent-based-model-development) and ran and posted a local test above from before and after development to ensure performance is monitored.

### Checks by code reviewer(s)

- [ ] Corresponding issue(s) is/are linked and addressed
- [ ] Code is clean of development artifacts (no deactivated or commented code lines, no debugging printouts, etc.)
- [ ] Appropriate **unit tests** have been added, CI passes, code coverage and performance is acceptable (did not decrease)
- [ ] No large data files added in the whole history of commits(files should in sum not exceed 100 KB, avoid PDFs, Word docs, etc.)
- [ ] On merge, add 2-5 lines with the changes (main added features, changed items, or corrected bugs) to the merge-commit-message. This can be taken from the **briefly-list-the-changes** above (best case) or the separate commit messages (worst case).
