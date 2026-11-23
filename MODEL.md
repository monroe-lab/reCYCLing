# reCYCLing model specification

This document specifies the forward-time model implemented by reCYCLing, in enough detail to reimplement it. It is the reference for methods sections that cite the package rather than restating the model.

Everything below is read out of the source as of version 1.9.2. Source locations are given in parentheses: `R/run_sim.R` holds the user-facing entry points and the pure-R backend; `src/sim_engine.cpp` holds the C++ engine. Where the source is ambiguous or where two code paths disagree, this document says so explicitly rather than choosing a reading.

The vignettes (`vignettes/reCYCLing-intro.Rmd`, `vignettes/equilibrium-dynamics.Rmd`, `vignettes/mechanisms-conversion-and-hotspots.Rmd`) are tutorials and are not duplicated here.

---

## 1. Data contract

A simulated genome is a `data.table` with one row per monomer and four columns (`R/run_sim.R`, `sim_genome()` return block):

| Column | Type | Meaning |
|---|---|---|
| `chrom` | integer | Chromosome index, `1..K`, gapless |
| `num` | integer | Position of the monomer within its chromosome, `1..n_c`, gapless and in array order |
| `bponly` | character | The monomer sequence, gapless, over the alphabet `A/C/G/T` |
| `dir` | character | Orientation label, `"+"` or `"-"` |

Rows are emitted chromosome by chromosome, in increasing `num` within each chromosome (`src/sim_engine.cpp`, output block of `sim_genome_cpp`).

The table carries four attributes (`R/run_sim.R`): `array_sizes` (integer vector of length `K`, final monomer count per chromosome), `total_gens` (generations actually run), `n_exchange` (count of completed inter-chromosomal exchange events), `hit_hard_cap` (logical).

`sim_array()` returns a list, not a bare table. Its `monomers` element is the same contract minus `chrom`: columns `num`, `bponly`, `dir` (`R/run_sim.R`, `.make_unit_dt()` and the C++ path's inline `data.table()` construction). The other elements are `pairs` (output of `pairs_identical()`, or an empty `data.table` when `compute_pairs = FALSE`), and the per-generation trajectories `L_vec`, `H_vec`, `N_vec`.

**Orientation is a label, not a transformation.** When a mechanism inserts an inverted chunk, the engine reverses the order of the monomers in the chunk and flips each monomer's `dir` sign; it does **not** reverse-complement `bponly` (`src/sim_engine.cpp`, the `invert` branches in `apply_within_array_moves`, `sim_core_cpp`, and the translocation block of `sim_genome_cpp`). Two monomers are therefore "identical" iff their stored `bponly` strings match, regardless of `dir`. Real data must be ingested under the same convention for `pairs_identical()` to mean the same thing on both.

**Real data enters the same contract** via `assemble_genome_dt()` (`R/assemble_genome.R`), which reads per-chromosome monomer FASTAs for one taxon and one haplotype, parses `chrom`, `num` and `dir` out of the FASTA headers by field position, and sets `bponly` to the header sequence with `-` characters stripped. Two differences from the simulated contract are worth stating: it returns a fifth column `aligned` (the gapped sequence, retained for divergence computation), and its `num` column is `numeric` rather than `integer`. Its `chrom` is the rank of the chromosome id within `sort(unique(chrom_id))`, which is a lexicographic sort of character ids — so unpadded ids (`Chr1`, `Chr2`, ... `Chr10`) will rank `Chr10` before `Chr2`.

---

## 2. Initialization

`sim_genome()` builds one ancestral monomer and seeds every chromosome from it (`R/run_sim.R`, then `src/sim_engine.cpp` seeding loop).

The ancestral monomer is `init_l` bases long (default 178). With `init_sequence_type = "random"` (default) each base is drawn i.i.d. uniformly from `c("A","C","G","T")` (`R/init_array.R`). With `init_sequence_type = "given"` the caller supplies `ancestor_seq`, which must have `nchar() == init_l`.

Chromosome `k` is seeded with exactly `target_sizes[k]` copies of that single shared monomer. All `sum(target_sizes)` monomers in the genome are identical at t = 0, and every chromosome starts exactly at its own target size. There is no growth phase: the genome is a "fully grown, zero-divergence" array that then ages in place for exactly `n_generations` generations of governed turnover.

Because all chromosomes descend from one monomer, cross-chromosome identical monomers can only arise from (i) inter-chromosomal exchange, (ii) inter-chromosomal gene conversion, or (iii) the vanishing chance of independent convergence — never from shared ancestry that persists, since divergence accumulates independently on each chromosome.

Monomer length is fixed for the whole run under the C++ engine: no mechanism inserts or deletes bases (section 4.7).

`sim_array()` seeds a single array with `init_k0` identical copies (default 10) of one ancestral monomer of length `init_l`, and then *grows* it, stopping when the array reaches `max_units` monomers or `max_generations` generations, whichever comes first (`R/run_sim.R`, `sim_array()` defaults; `src/sim_engine.cpp`, the `while` condition in `sim_core_cpp`). This is a different experimental setup from `sim_genome()` — grow-to-size versus age-at-size — using the same mechanism set.

---

## 3. Per-generation update order

Order matters; this is the exact sequence.

### 3.1 `sim_genome()` (`src/sim_engine.cpp`, main loop of `sim_genome_cpp`)

For each of `n_generations` generations:

1. **Guards.** For each chromosome: set the `hit_hard_cap` flag if its size exceeds `hard_caps[c]` (informational only; the run continues), and if it has reached size 0, re-seed it with `RESEED_K0 = 10` ancestral monomers in a new lineage family.
2. **Global governor.** Compute one genome-wide per-monomer deletion probability `p_del_global` from the total monomer count at the start of this generation (section 4.6).
3. **Within-chromosome moves**, applied to chromosomes in index order `c = 1..K`, via the shared helper `apply_within_array_moves()`. Within each chromosome the sub-order is fixed:
   1. local (tandem) duplication;
   2. long-range (distal) within-chromosome duplication;
   3. chunk deletion (using this chromosome's governed rate `p_del_c`);
   4. gene conversion, within-chromosome path;
   5. terminal deletion — **not reachable from `sim_genome()`**, which hard-codes `p_distal_del = 0.0`.
4. **Inter-chromosomal exchange**, at most one event for the whole genome.
5. **Inter-chromosomal gene conversion**, scanning every monomer of every chromosome.
6. **Clock tick.** Increment every monomer's pending-mutation counter by one; if `hotspot_alpha > 0`, also roll each monomer's lineage-family break (section 4.4).

After the loop, all pending mutations are materialized on every monomer and the output table is built.

### 3.2 `sim_array()`, C++ backend (`src/sim_engine.cpp`, `sim_core_cpp`)

Per generation, on the single array: hard-cap check and extinction re-seed → local duplication → distal duplication → chunk deletion → gene conversion → terminal deletion → increment pending-mutation counters. Mutations are materialized at the end of the run.

### 3.3 `sim_array()`, R backend (`R/run_sim.R`, `.step_generation()`)

Per generation: local duplication → distal duplication → chunk deletion → per-base mutation, applied immediately. Gene conversion and terminal deletion are not implemented in the R backend; requesting either with `backend = "r"` raises an error (`.run_one_replicate()`).

---

## 4. Mechanisms

### 4.0 Chunk-size distributions

Every chunk size, donor separation and conversion tract length is drawn from a distribution specified as a list with a `type` field (`R/run_sim.R`, `.sample_chunk_size()`; `src/sim_engine.cpp`, `sample_chunk_size()`), one of:

`list(type="fixed", value=)`, `list(type="poisson", lambda=)`, `list(type="normal", mean=, sd=)`, `list(type="geom", prob=)`, `list(type="unif", min=, max=)`, `list(type="gamma", shape=, scale=)` (or `rate=` instead of `scale=`).

Every draw is rounded to an integer, floored at 1, and capped at a mechanism-specific maximum `max_k`. Two consequences follow and apply everywhere below. First, a nominal `gamma(shape=1, scale=20)` donor separation is *not* exactly Exponential(mean 20): it is that distribution rounded and floored at 1, which moves a small amount of mass (`P(X < 0.5) ≈ 0.025`) up to 1. Second, all chunks extend *rightward* from their sampled start position and are truncated at the end of the array, so realized chunk sizes are shorter than nominal near the array's right end.

### 4.1 Local (tandem) duplication

Rate parameter `p_local_dup`, **per monomer per generation** (default `1.5e-4` in both `sim_array()` and `sim_genome()`).

Each monomer independently triggers with probability `p_local_dup`. For a trigger at position `i`, a chunk size is drawn from `local_dist` (default `gamma(shape=2, scale=15)`, mean 30 monomers), capped at the number of monomers from `i` to the end of the array. The chunk `[i, i+size-1]` is copied and the copy is inserted **immediately after the source chunk**, preserving monomer order and every `dir`. There is no inversion on local duplication.

Triggers are collected first against the array as it stands at the start of the phase, then applied in **reverse index order**, so insertions do not shift the indices of pending triggers (`src/sim_engine.cpp`, the reverse loop over `triggered`; the R backend does the same with `rev(idxs)`).

### 4.2 Long-range (distal) within-chromosome duplication

Rate parameter `p_distal_dup`, **per chromosome per generation** in `sim_genome()` and **per array per generation** in `sim_array()` (default 0 in both, i.e. off). It is the same code in both cases; the difference is only that `sim_genome()` runs it once per chromosome per generation. At most one event per chromosome per generation.

On an event: a start position is drawn uniformly over the array; a chunk size is drawn from `distal_dist` (default `gamma(shape=2, scale=500)`, mean 1000 monomers) and capped at the distance to the array end; the chunk is copied, inverted with probability `p_invert_distal`, and inserted at a position drawn uniformly over the `k+1` insertion slots of the array **as sized before the insertion** — that is, the target position is uniform over the whole chromosome, not restricted to be far from the source.

### 4.3 Inter-chromosomal exchange

Rate parameter `p_exchange`, **per genome per generation** (default 0, i.e. off). `sim_genome()` only; it does not exist in `sim_array()`. At most one event per generation, and it requires `K >= 2`.

The R argument `p_exchange` is passed to the C++ argument `p_translocation`, and `exchange_dist` to `transloc_dist` (`R/run_sim.R`, the `sim_genome_cpp()` call). The two names denote the same mechanism.

On an event: a donor chromosome is drawn uniformly over `1..K`; a recipient is drawn uniformly over `1..K` and, if it equals the donor, replaced by `(donor + 1) mod K`. **This makes the recipient non-uniform**: conditional on the donor, the chromosome at index `donor+1 (mod K)` is chosen with probability `2/K` and each other non-donor chromosome with probability `1/K` (`src/sim_engine.cpp`, `if (recip == donor) recip = (recip + 1) % K;`). This is a documented property of the code, not a stated modelling intent.

If the recipient is already at or above its hard cap the event is abandoned (the donor may still donate to other recipients on later generations); otherwise a start position is drawn uniformly on the donor, a chunk size is drawn from `exchange_dist`, capped at the distance to the donor's end, and the chunk is **copied** (not moved) into the recipient at a uniformly drawn insertion slot, inverted with probability `p_invert_exchange`. The donor keeps its copy.

Exchange is therefore modelled as duplicative, conversion-like transfer. In single-genome data reciprocal and non-reciprocal exchange are not distinguishable, and the model does not attempt to.

`exchange_dist` defaults to `NULL`, in which case `sim_genome()` substitutes `gamma(shape=1, scale=exchange_chunk_mean)` with `exchange_chunk_mean = 500` (`R/run_sim.R`) — an exponential-shaped block-size distribution with nominal mean 500 monomers, in contrast to the `shape=2` distributions used for local and distal duplication.

Completed events are counted and returned as `attr(g, "n_exchange")`.

### 4.4 Duplication-hotspot feedback

Rate parameter `hotspot_alpha >= 0` (default 0, off). `sim_genome()` only. It modifies where local duplications land, not how many there are.

When `hotspot_alpha > 0`, the per-monomer local-duplication trigger probability becomes

```
w_i = (size of monomer i's lineage family on this chromosome) ^ alpha
p_i = p_local_dup * k * w_i / sum_j(w_j)
```

where `k` is the current chromosome size (`src/sim_engine.cpp`, the `hotspot_alpha > 0.0` branch of `apply_within_array_moves`). Family weights are recomputed from scratch every generation, for each chromosome separately.

Three properties follow directly from that normalization:

- **Mean-preserving.** `sum_i p_i = p_local_dup * k` for any `alpha`, so the expected number of local-duplication triggers per chromosome per generation is `p_local_dup * k` — unchanged by `alpha`. `alpha` redistributes duplication events toward large families without adding duplication mass, so it is not confounded with `p_local_dup` and the size governor's balance arithmetic (section 4.6) is untouched. The one exception: `p_i` is used in a Bernoulli draw, so if `p_i > 1` it saturates at 1 and the realized mean falls below `p_local_dup * k`. This only occurs at extreme `alpha`.
- **`alpha = 0` recovers the uniform model exactly.** The code branches on `hotspot_alpha > 0.0` and takes the plain uniform path otherwise, consuming exactly `k` uniform deviates either way, and the family-break roll in the clock tick is skipped entirely when `alpha == 0`. Runs with `alpha = 0` are therefore bit-identical to runs of the engine with the mechanism absent, under the same seed.
- **Rich-get-richer.** Larger families get more duplications, which make them larger.

**Correction on "identical family".** `w_i` is not computed over exact-sequence identity classes. It is computed over a **lineage family id** carried on each monomer (`struct Monomer { ... long long fam; }`). Copies inherit their source's id — local, distal, exchanged and gene-conversion-recipient copies all join their source's or donor's family. Each monomer independently leaves its family, receiving a fresh singleton id, with per-generation probability

```
p_fam_break = 1 - (1 - mu_total)^L
```

i.e. the exact probability that at least one substitution lands anywhere in a monomer of length `L` in one generation (`src/sim_engine.cpp`, `p_fam_break` definition and the clock-tick loop). The source comment gives the reason: under deferred mutation the stored sequences are stale, so hashing them would see the seeded ancestral family as identical forever. Family membership therefore tracks exact-sequence identity **in distribution**, at the correct rate, but not in per-base realization: a monomer can break its family while its realized sequence is unchanged, and vice versa. Seeding assigns one family per chromosome (`fam = c`), so families are within-chromosome by construction from t = 0.

### 4.5 Gene conversion

Rate parameter `p_conversion`, **per monomer per generation** (default 0, off). Available in `sim_genome()` and in `sim_array()`'s C++ backend only; the R backend errors out if it is requested.

A conversion event **overwrites a recipient monomer's sequence with a donor's sequence**. Monomer count does not change, chromosome sizes do not change, and no monomer moves. The mechanism homogenizes without adding units. The recipient's `dir` is **not** changed (orientation is treated as a property of the locus, not of the converted tract). The recipient's lineage family id is set to the donor's.

In `sim_genome()` the single rate is split into two independent per-monomer Bernoulli paths (`src/sim_engine.cpp`, `p_conv_within` / `p_conv_xchrom`):

```
p_conv_within = p_conversion * (1 - p_conv_interchrom)
p_conv_xchrom = p_conversion * p_conv_interchrom
```

`p_conv_interchrom` defaults to 0 (all conversion within-chromosome) and must lie in `[0, 1]`. The split is on the *rate*, implemented as two separate draws, so `p_conv_interchrom` is the expected fraction of conversion events that are inter-chromosomal, not a literal partition of a single donor draw.

**Within-chromosome path** (inside `apply_within_array_moves`): each monomer `i` is a recipient with probability `p_conv_within`. A separation `off` is drawn from `conv_dist`, default `gamma(shape=1, scale=20)` — an exponential with nominal mean 20 monomers, subject to the rounding and flooring of section 4.0 and capped at `k-1`. A direction sign is drawn `+1` or `-1` with probability 1/2 each, giving donor index `i + sign*off`. If that falls off either end it is **reflected** to the other side (`i + off` or `i - off`); if it is still out of range, or equals `i`, the event is abandoned. A tract length is drawn from `conv_tract`, default `gamma(shape=1, scale=1)` — an exponential with nominal mean 1 monomer, rounded and floored at 1, and capped at the room remaining to the right of both recipient and donor. The tract is copied monomer by monomer, donor position `donor+j` onto recipient position `i+j`.

**Inter-chromosomal path** (genome level, step 5): each monomer of each chromosome is a recipient with probability `p_conv_xchrom`. The donor chromosome is drawn uniformly from the `K-1` chromosomes other than the recipient's, and the donor position is drawn **uniformly over that chromosome** — `conv_dist` is not used on this path, since separation along an array is undefined across chromosomes. Tract length is drawn from `conv_tract` as above. After an event the scan index skips the converted tract (`i += tract - 1`), so within one generation the inter-chromosomal path cannot convert the same tract twice; the within-chromosome path has no such skip.

Two side effects of the implementation that a reimplementation must match. Before the sequence is copied, the *donor* is materialized to its current sequence (section 4.7); the *recipient's* own pending mutations are discarded, not applied — `arr[rpos].pending_gens = 0` — which is correct, since the recipient's sequence is being overwritten anyway. And conversion does not touch `num`, so the array's coordinate structure is unchanged.

### 4.6 Deletion and the genome-wide size governor

There are two deletion mechanisms in the engine and `sim_genome()` exposes only the governed one.

**Chunk deletion** is triggered **per monomer per generation**. For a trigger at `i`, a chunk size is drawn from `del_dist` (default `gamma(shape=2, scale=15)`, mean 30 monomers) and capped at the distance to the array end, and monomers `[i, i+size-1]` are erased. As with duplication, triggers are collected first and applied in reverse index order. In `sim_array()` the rate is the user parameter `p_del_chunk` (default 0). In `sim_genome()` there is no user-facing deletion rate at all: the rate is set by the governor every generation.

**Terminal deletion** (`p_terminal_del` in `sim_array()`, `p_distal_del` in the C++ signature) is **per array per generation**, at most one event, and removes a chunk of size drawn from `distal_dist` from the **end** of the array. It is the mirror of distal duplication and models telomere-proximal truncation. It is off by default, is C++-only, and is **unreachable from `sim_genome()`**, which passes `p_distal_del = 0.0` unconditionally.

**The governor** (`src/sim_engine.cpp`, the block above the per-chromosome loop in `sim_genome_cpp`). Each generation, one genome-wide per-monomer deletion probability is computed:

```
p_del_balance = p_local_dup * mean_local / mean_del
              + p_distal_dup * mean_distal * K / (total_sz * mean_del)

err_raw = (total_sz - total_target) / total_target
err     = err_raw - band   if err_raw >  band
        = err_raw + band   if err_raw < -band
        = 0                otherwise

p_del_global = clamp(p_del_balance + gov_gain * err / mean_del, 0, 0.5)
```

where `total_sz` is the genome-wide monomer count at the start of the generation, `total_target = sum(target_sizes)`, `band = size_band` (default 0.10), `gov_gain = 0.5` (a hard-coded constant in the C++ source), and `mean_local`, `mean_distal`, `mean_del` are the analytic means of the respective chunk-size distributions computed by `dist_mean()`.

Then, for each chromosome, a weak per-chromosome restoring term is added on top, using the chromosome's own size and target and the same deadband shape:

```
e_raw  = (sz_c - target_c) / target_c, banded the same way to give e
p_del_c = clamp(p_del_global + perchrom_gain * e / mean_del, 0, 0.5)
```

with `perchrom_gain = 0.01`, fifty times weaker than `gov_gain`, also hard-coded. Because the per-chromosome targets sum to the total target, the two terms pull in the same direction. The intent recorded in the source is that a chromosome keeps its characteristic size over the long run while a transient swell from an exchange event erodes slowly, at the background turnover rate, rather than being cancelled immediately.

Finally, a chromosome at or above its hard cap has `p_local_dup` and `p_distal_dup` set to 0 for that generation — it stops growing but keeps deleting — while the rest of the genome carries the size. `hard_caps` defaults to `ceiling(3 * max(target_sizes))`, applied to every chromosome (`R/run_sim.R`).

Four points where a plain-language description of the governor can go wrong:

- The balance term matches **expected** duplication mass, computed from the distribution **means** and the nominal rates. It does not observe the duplications that actually happened in the generation. Realized net drift per generation is therefore nonzero and the deadband absorbs it.
- `size_band` is a **deadband**, not a tolerance. Inside `±size_band` of the total target there is *no* proportional correction at all and the total floats freely on the net-zero balance. The proportional pull only switches on once the total has already drifted outside the band. The code does not guarantee that the total stays within `±10%`.
- The balance term accounts for local duplication and within-chromosome distal duplication only. **Inter-chromosomal exchange mass is not in it**, even though exchange is duplicative and adds monomers to the genome. With `p_exchange > 0` there is an uncompensated source of growth that only the band correction removes. Gene conversion needs no term, since it does not change monomer count.
- Both `p_del_global` and `p_del_c` are clamped into `[0, 0.5]`.

### 4.7 The per-base substitution clock

Rate parameter `mu_total`, **per base per generation** (default `1e-4` in `sim_array()`, `5e-5` in `sim_genome()`).

The substitution model is Jukes-Cantor: with probability `mu_total` per base per generation, a base changes to one of the three other bases chosen uniformly.

The C++ engine does not apply this base by base every generation. Each monomer carries an integer counter `pending_gens` of mutation rounds not yet applied (`struct Monomer`). The clock tick at the end of every generation increments the counter on every monomer. `materialize()` then collapses `g` deferred generations in one pass using the Jukes-Cantor batch transition probabilities:

```
r      = (1 - 4*mu_total/3)^g
p_same = 1/4 + 3/4 * r
```

Each base is redrawn to a uniformly chosen different base with probability `1 - p_same`, and the counter is reset to 0. For `g = 1` this gives `p_same = 1 - mu_total`, matching the per-generation definition exactly; for `g > 1` it correctly accounts for multiple hits and back-mutation. The scheme is exactly equivalent **in distribution** to applying one generation at a time; it is not equivalent realization-by-realization under a given seed.

Materialization is lazy and is forced only where a sequence is about to be read or copied:

- on the source monomers of a local duplication, before the copy is made;
- on the source monomers of a distal duplication, before the copy;
- on the donor monomers of an inter-chromosomal exchange, before the copy;
- on the donor monomer of each gene-conversion step, before the overwrite;
- on every monomer at the end of the run, before the output table is built.

The consequence that matters for interpretation: at the instant a chunk is duplicated, source and copy are exactly identical, and from then on they accumulate substitutions independently. Divergence between two monomers is a function of the time since their most recent copy event, which is what the model is for.

What the substitution clock does **not** model: insertions or deletions of bases (so monomer length is constant under the C++ engine — the R backend's `p_ins` and `p_del_bp` are the only path to indels, and neither `sim_genome()` nor the C++ `sim_array()` backend supports them); among-site rate variation; transition/transversion bias or any departure from Jukes-Cantor; base-composition bias (the ancestral monomer is uniform-random unless supplied); selection of any kind; and any population-genetic process — there is no population, no diploidy, no reproduction, no recombination between haplotypes, and no genetic drift on allele frequencies. A reCYCLing run is one haplotype's array content evolving by mutation and rearrangement.

---

## 5. Time units

A generation in this model has no organismal meaning. The code contains no reproduction, no population size, and no life-history process; a "generation" is simply one pass of the update loop in section 3.

The only quantity that gives the generation an absolute scale is the substitution clock. With `mu_total` fixed, divergence from the array consensus grows as `mu_total × effective turnover age`, so `n_generations` is the divergence dial, and any conversion between simulated generations and calendar time has to be made externally by fixing `mu_total` to a per-base per-generation substitution rate on an independent basis. All other rates (`p_local_dup`, `p_distal_dup`, `p_exchange`, `p_conversion`) are therefore only interpretable *relative to* `mu_total`.

---

## 6. Parameter glossary

Defaults are those in the exported R function signatures (`R/run_sim.R`) unless noted.

### `sim_genome()`

| Parameter | Controls | Units / per what | Default |
|---|---|---|---|
| `K` | Number of chromosome arrays | count | required |
| `target_sizes` | Seed size and governor target of each chromosome | monomers, length `K` | required |
| `n_generations` | Length of the run; the divergence dial | generations | required |
| `mu_total` | Substitution clock | per base per generation | `5e-5` |
| `p_local_dup` | Local (tandem) duplication trigger | per monomer per generation | `1.5e-4` |
| `local_dist` | Local duplication chunk size | monomers | `gamma(shape=2, scale=15)` |
| `p_distal_dup` | Long-range within-chromosome duplication | per chromosome per generation | `0` |
| `distal_dist` | Distal duplication chunk size | monomers | `gamma(shape=2, scale=500)` |
| `p_exchange` | Inter-chromosomal exchange | per genome per generation | `0` |
| `exchange_dist` | Exchanged block size | monomers | `NULL` → `gamma(shape=1, scale=exchange_chunk_mean)` |
| `exchange_chunk_mean` | Mean exchanged block size when `exchange_dist` is `NULL` | monomers | `500` |
| `p_invert_distal` | Inversion of a distal-duplicated chunk | probability per event | `0.22` |
| `p_invert_exchange` | Inversion of an exchanged chunk | probability per event | `0.35` |
| `del_dist` | Governor deletion chunk size | monomers | `gamma(shape=2, scale=15)` |
| `p_conversion` | Gene conversion (recipient trigger) | per monomer per generation | `0` |
| `conv_dist` | Within-chromosome donor–recipient separation | monomers | `gamma(shape=1, scale=20)` |
| `conv_tract` | Conversion tract length | monomers | `gamma(shape=1, scale=1)` |
| `p_conv_interchrom` | Fraction of the conversion rate routed to an inter-chromosomal donor | fraction in `[0,1]` | `0` |
| `hotspot_alpha` | Exponent of the family-size duplication weight | dimensionless, `>= 0` | `0` |
| `size_band` | Governor deadband, as a fraction of total target | fraction in `(0,1)` | `0.10` |
| `hard_caps` | Per-chromosome memory guard | monomers | `NULL` → `ceiling(3 * max(target_sizes))` |
| `init_l` | Monomer length | bp | `178` |
| `init_sequence_type` | Ancestral monomer source | `"random"` / `"given"` | `"random"` |
| `ancestor_seq` | Supplied ancestral monomer | character | `NULL` |
| `verbose` | Print chromosome sizes every 500 generations | logical | `FALSE` |

### `sim_array()` (additional or differing)

| Parameter | Controls | Units / per what | Default |
|---|---|---|---|
| `max_units` | Stop when the array reaches this size | monomers | `20000` |
| `max_generations` | Generation cap | generations | `Inf` |
| `hard_cap` | Safety stop | monomers | `40000` |
| `init_k0` | Initial copy number | monomers | `10` |
| `mu_total` | Substitution clock | per base per generation | `1e-4` |
| `p_local_dup` | Local duplication trigger | per monomer per generation | `1.5e-4` |
| `p_distal_dup` | Distal duplication | per **array** per generation | `0` |
| `p_invert_distal` | Inversion of a distal-duplicated chunk | probability per event | `0.22` |
| `p_del_chunk` | Chunk deletion trigger (ungoverned) | per monomer per generation | `0` |
| `p_terminal_del` | Terminal truncation from the array end, size `~ distal_dist` | per array per generation, C++ only | `0` |
| `p_conversion` | Gene conversion | per monomer per generation, C++ only | `0` |
| `p_sub`, `p_ins`, `p_del_bp` | Relative weights of substitution / base insertion / base deletion | relative weights, R backend only | `1`, `0`, `0` |
| `backend` | Engine selection | `"auto"` / `"cpp"` / `"r"` | `"auto"` |
| `compute_pairs` | Whether to enumerate identical pairs (O(n²)) | logical | `TRUE` |

### Hard-coded constants (`src/sim_engine.cpp`)

| Constant | Value | Role |
|---|---|---|
| `gov_gain` | `0.5` | Proportional gain of the global size governor |
| `perchrom_gain` | `0.01` | Proportional gain of the weak per-chromosome restoring term |
| deletion-probability clamp | `[0, 0.5]` | Applied to both `p_del_global` and `p_del_c` |
| `RESEED_K0` | `10` | Monomers used to re-seed a chromosome that deletes to extinction |

`p_invert_distal = 0.22` and `p_invert_exchange = 0.35` are **defaults in the exported R signatures** of both `sim_array()` (for `p_invert_distal`) and `sim_genome()`, not values that a caller must supply. A caller may override them. The internal R helper `.distal_duplication()` carries its own fallback default of `0.5`, but `sim_array()` always passes an explicit value, so that fallback is never used in practice. Note also that no inversion is applied by local duplication or by gene conversion.

---

## 7. The two backends

`sim_genome()` is C++-only. `sim_array()` has both, selected by `backend`: `"auto"` (default) uses C++ unless base-level indels are requested (`p_ins > 0` or `p_del_bp > 0`); `"cpp"` forces C++ and errors if indels are requested; `"r"` forces the R backend.

| | C++ (`src/sim_engine.cpp`) | R (`R/run_sim.R`) |
|---|---|---|
| Speed | ~100× faster (stated in `README.md`; not measured here) | baseline |
| Mutation | Deferred; Jukes-Cantor batch materialization | Applied base by base every generation |
| Substitutions | yes | yes |
| Base insertions / deletions (`p_ins`, `p_del_bp`) | **no** | yes |
| Gene conversion (`p_conversion`) | yes | **no** — errors |
| Terminal deletion (`p_terminal_del`) | yes | **no** — errors |
| Hotspot feedback (`hotspot_alpha`) | yes, `sim_genome()` only | n/a |
| `H_vec`, `N_vec` | final value only; all earlier entries `NA_real_` | full per-generation trajectory |
| Sequence storage | integer-encoded bases (0=A,1=C,2=G,3=T) | character strings |

The two backends are **not** seed-compatible. They consume random deviates in different orders and quantities, and deferred versus per-generation mutation differs realization-by-realization even though it agrees in distribution.

One further difference is a genuine off-by-one that the source does not reconcile, and that a reimplementer should be aware of. The C++ loop is `t = 0; while (t < max_t) { t++; ...apply moves... }`, which performs `max_t` update steps. The R loop is `t = 1; while (... && t < max_generations) { t <- t + 1; ...apply moves... }`, which performs `max_generations - 1` update steps. For the same `max_generations = N`, the C++ backend ages the array `N` generations and the R backend `N - 1`. The `L_vec` lengths differ correspondingly. Nothing in the source states which is intended.

The extinction handlers also differ slightly: on reaching size 0, the C++ engine re-seeds and `continue`s, skipping that generation's moves, whereas the R backend re-seeds and then applies that generation's moves normally. Both reset the generation counter, both warn. Under `sim_genome()`'s governor, extinction is effectively unreachable.
