# To-Do List

## FragmentRetro

- [x] Finish `Retrosynthesis` class
- [x] Create a `Solution` class to visualize retrosynthesis "solutions".
- [x] Rewrite `is_substructure_BBs` such that BBs change at every stage (take out unnecessary BBs).
- [x] Write tests for the `CompoundFilter` class (with or without filter should work the same for substructure matcher)
- [x] Write tests for the `RetrosynthesisSolution` class
- [x] Write tests for the parallelization of substructure matcher
- [x] Write tests for `BRICSFragmenter`.
- [x] Write tests for `get_combination_smiles`.
- [ ] Ignore chirality at the neighbor of dummy atoms for `is_strict_substructure` (use both @ and @@).
- [ ] Confirm that `SubstructureMatcher.addH_to_wildcard_neighbors` is the most efficient way to do `is_strict_substructure`.
- [ ] Write tests for the `Retrosynthesis` class

## New Fragmentation Rules

- [ ] Come up with ring fragmentation rules that covers enough common types of rings in drug-like compounds.
