# Releasing paper-ARSA

A release is a published GitHub release with a tag `vX.Y.Z`. Zenodo archives it, gives it a DOI, and takes the record's metadata (title, author, description, keywords, license, version) from `CITATION.cff`. The paper cites that DOI. Make the first release, `v1.0.0`, when the paper is final; a revision of the paper that changes a figure or a number gets a new release.

## Before the first release

1. **Make the repository public.** Zenodo archives only public repositories:

   ```bash
   gh repo edit chongchonghe/paper-ARSA --visibility public --accept-visibility-change-consequences
   ```

2. **Enable the repository on Zenodo.** At <https://zenodo.org/account/settings/github/> (already linked to GitHub for Honghu.jl), find `chongchonghe/paper-ARSA` and turn its switch on. Zenodo archives only releases published after this.
3. **Decide the authors of the record.** `CITATION.cff` lists the author of the code; add the co-authors of the paper if they should be on the record too.

## Each release

1. **Rebuild everything** from a clean clone of the commit to be released, so that every header names a commit of this repository and none says `-dirty`:

   ```bash
   git clone https://github.com/chongchonghe/paper-ARSA && cd paper-ARSA
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ./run_all.sh
   git status        # only the headers of the printouts and the run times should change
   ```

   Compare the regenerated data with those in git; commit them if they agree, and check that the paper quotes the same numbers. Then run `./copy_to_paper.sh` and check that the paper compiles with the figures.
2. **Update `CITATION.cff`**: add `version: X.Y.Z` and `date-released: YYYY-MM-DD`, and the paper's DOI or arXiv identifier in `preferred-citation` once it has one. Check the file:

   ```bash
   uvx --from cffconvert cffconvert --validate -i CITATION.cff   # the metadata Zenodo will read
   uvx --from cffconvert cffconvert -f zenodo -i CITATION.cff    # ... and what Zenodo will make of it
   ```

3. **Commit, push and publish the release**:

   ```bash
   v=1.0.0
   gh release create v$v --target main --title "paper-ARSA v$v" \
       --notes "Scripts and data of He, Krumholz & Wibking, \"An Adaptive Reduced Speed of Light Approximation\"."
   ```

   Do not mark it as a draft: Zenodo archives published releases.
4. **Check the record** at <https://zenodo.org/account/settings/github/repository/chongchonghe/paper-ARSA>. Edit it on Zenodo to add the paper as a related work ("is supplement to", with its DOI or arXiv identifier); `CITATION.cff` cannot pass that link to Zenodo, and the record can be edited without a new version.
5. **Cite the DOI.** Put the concept DOI (the one that resolves to the latest version) in the paper's code and data availability statement and in Sec. 5 next to Honghu.jl, as `doi` in `CITATION.cff`, and as a badge at the top of `README.md`. In the outer repository, commit the new submodule commit.
