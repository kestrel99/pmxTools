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

We checked all 3 reverse dependencies (PKNCA, mrgsolve and rxode2; all list
pmxTools in Suggests) against pmxTools 1.6 and found no new problems:

* PKNCA 0.12.1 and mrgsolve 2.0.1: R CMD check, no problems related to
  pmxTools.
* rxode2 5.1.7.1: its tests that compare against pmxTools' closed-form
  solutions (test-solComp.R, run with NOT_CRAN=true) all pass.
