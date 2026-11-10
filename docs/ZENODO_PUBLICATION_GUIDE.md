# GitHub and Zenodo release sequence

## GitHub repository

Target repository: `https://github.com/Ainiyar/SAC-SQ-Plus`.

1. Create the repository as **Public**.
2. Push the complete contents of this release candidate.
3. Confirm the root contains `LICENSE`, `DATA_LICENSE.md`, `CITATION.cff`, and `.zenodo.json`.
4. Create a GitHub Release tagged `v1.0.0` with the included release notes.

## Zenodo

1. Sign in to Zenodo using the GitHub account.
2. Authorize Zenodo to access public repositories.
3. Enable `Ainiyar/SAC-SQ-Plus` in the Zenodo GitHub integration.
4. Create the GitHub Release `v1.0.0` only after the repository is enabled. Zenodo then archives that release and assigns a version DOI; the record also provides a concept DOI for the software across versions.
5. Open the resulting Zenodo record and make sure the mixed-license declaration lists **MIT** for source code/scripts and **CC BY 4.0** for the authors' derived data, frozen cost boxes, processed windows, experiment tables, and figures. Zenodo supports declaring multiple licenses for mixed-license uploads; `.zenodo.json` provides the MIT software license, while `DATA_LICENSE.md` documents the CC BY 4.0 data scope. Add the second license in the Zenodo record's license interface if it was not imported automatically.
6. Add the Zenodo version DOI to the manuscript's Data and Code Availability statement and the repository's citation metadata. The concept DOI may also be cited when referring to the software family independent of a specific release.

Official instructions:
- GitHub: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content
- Zenodo GitHub archive: https://help.zenodo.org/docs/github/archive-software/
- Zenodo mixed licenses: https://help.zenodo.org/docs/deposit/describe-records/licenses/
- Zenodo `.zenodo.json`: https://help.zenodo.org/docs/github/describe-software/zenodo-json/

## DOI placement in R400

After the DOI is minted, use the version-specific DOI in the manuscript's Code Availability / Data Availability statement and cite the GitHub repository as the live development location. Do not invent or pre-fill a DOI before Zenodo creates the record.
