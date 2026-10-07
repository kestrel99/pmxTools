## Submission

This is an update from pmxTools 1.5 to 1.6. It adds `cut_quantile()` and
`datamap()`, changes the `dgr_table()` interface (the old form still works),
and fixes `read_nm()` and a load-time import warning. Full details are in
NEWS.md.

The reference URL for Bertrand & Mentre (2008) now points to an
Internet Archive snapshot, because the original host's TLS certificate has
expired.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local Windows 11, R 4.6.1
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), macos-latest
  (release), windows-latest (release)

## Reverse dependencies

pmxTools has three reverse dependencies, all in Suggests: PKNCA, mrgsolve and
rxode2.

TODO: record revdep check results before submitting.
