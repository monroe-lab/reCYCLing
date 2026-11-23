#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <cstdint>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Monomer representation: integer-encoded bases + deferred mutation counter
// ---------------------------------------------------------------------------
// Bases: 0=A, 1=C, 2=G, 3=T
// pending_gens: number of mutation rounds not yet applied
// dir: 1 = "+", -1 = "-"

struct Monomer {
    std::vector<int> seq;
    int pending_gens;
    int dir;
    // Lineage family id for the hotspot-feedback mechanism (g9). Copies inherit
    // their source's id; a monomer leaves its family (fresh singleton id) with
    // the per-generation probability that >=1 substitution lands anywhere in
    // the monomer -- exactly the rate at which true exact-sequence identity
    // breaks. Tracked this way because stored seqs are STALE under deferred
    // mutation (hashing them would see the seeded ancestor family as identical
    // forever). Maintained only when hotspot_alpha > 0; ignored otherwise.
    long long fam;

    Monomer() : pending_gens(0), dir(1), fam(0) {}
    Monomer(const std::vector<int>& s, int d) : seq(s), pending_gens(0), dir(d),
                                                fam(0) {}
};

// ---------------------------------------------------------------------------
// Chunk size sampling (mirrors R's .sample_chunk_size)
// ---------------------------------------------------------------------------

// Expected (mean) chunk size for a chunk-size distribution spec. Used by the
// genome steady-state governor to balance duplication mass against deletion
// mass (net-zero turnover). Falls back to 1.0 if the mean is undefined/<1.
double dist_mean(const List& dist) {
    std::string type = as<std::string>(dist["type"]);
    double m;
    if (type == "fixed") {
        m = as<double>(dist["value"]);
    } else if (type == "poisson") {
        m = as<double>(dist["lambda"]);
    } else if (type == "normal") {
        m = as<double>(dist["mean"]);
    } else if (type == "geom") {
        double p = as<double>(dist["prob"]);
        m = (p > 0.0) ? (1.0 - p) / p : 1.0;   // R::rgeom mean = (1-p)/p
    } else if (type == "unif") {
        m = 0.5 * (as<double>(dist["min"]) + as<double>(dist["max"]));
    } else if (type == "gamma") {
        double shape = as<double>(dist["shape"]);
        double scale_val = dist.containsElementNamed("scale")
            ? as<double>(dist["scale"]) : 1.0 / as<double>(dist["rate"]);
        m = shape * scale_val;
    } else {
        m = 1.0;
    }
    if (!(m >= 1.0)) m = 1.0;   // also catches NaN
    return m;
}

int sample_chunk_size(const List& dist, int max_k) {
    std::string type = as<std::string>(dist["type"]);
    int s;

    if (type == "fixed") {
        s = (int)as<double>(dist["value"]);
    } else if (type == "poisson") {
        s = (int)R::rpois(as<double>(dist["lambda"]));
    } else if (type == "normal") {
        s = (int)std::round(R::rnorm(as<double>(dist["mean"]),
                                      as<double>(dist["sd"])));
    } else if (type == "geom") {
        s = (int)R::rgeom(as<double>(dist["prob"]));
    } else if (type == "unif") {
        s = (int)std::round(R::runif(as<double>(dist["min"]),
                                      as<double>(dist["max"])));
    } else if (type == "gamma") {
        double shape = as<double>(dist["shape"]);
        double scale_val;
        if (dist.containsElementNamed("scale")) {
            scale_val = as<double>(dist["scale"]);
        } else {
            scale_val = 1.0 / as<double>(dist["rate"]);
        }
        s = (int)std::round(R::rgamma(shape, scale_val));
    } else {
        stop("Unknown chunk size distribution type: %s", type.c_str());
    }

    if (s < 1) s = 1;
    if (s > max_k) s = max_k;
    return s;
}

// ---------------------------------------------------------------------------
// Materialize pending mutations using Jukes-Cantor batch formula
// ---------------------------------------------------------------------------
// After g generations of JC substitution with per-base rate mu:
//   P(base unchanged) = 1/4 + 3/4 * (1 - 4*mu/3)^g
//   P(base -> specific other) = (1 - (1 - 4*mu/3)^g) / 4
//
// This is exactly equivalent in distribution to applying mutations
// one generation at a time.

void materialize(Monomer& m, double mu_total) {
    int g = m.pending_gens;
    if (g == 0 || mu_total <= 0.0) return;

    double r = std::pow(1.0 - 4.0 * mu_total / 3.0, (double)g);
    double p_same = 0.25 + 0.75 * r;

    int L = (int)m.seq.size();
    for (int i = 0; i < L; i++) {
        if (R::runif(0.0, 1.0) >= p_same) {
            // Change to a uniformly random different base
            int offset = 1 + (int)(R::runif(0.0, 1.0) * 3);
            if (offset > 3) offset = 3; // guard against runif == 1.0
            m.seq[i] = (m.seq[i] + offset) % 4;
        }
    }
    m.pending_gens = 0;
}

// ---------------------------------------------------------------------------
// Per-array within-chromosome moves (shared by single-array and multi-chrom).
// Operates IN PLACE on one array. Mirrors exactly the move logic inside
// sim_core_cpp so both engines behave identically per chromosome.
// (sim_core_cpp keeps its own inline copy unchanged for safety; this helper is
// used only by the new multi-chromosome engine sim_genome_cpp.)
// ---------------------------------------------------------------------------

void apply_within_array_moves(std::vector<Monomer>& arr,
                              double p_local_dup,
                              const List& local_dist,
                              double p_distal_dup,
                              const List& distal_dist,
                              double p_invert_distal,
                              double p_del_chunk,
                              const List& del_dist,
                              double p_distal_del,
                              double p_conversion,
                              const List& conv_dist,
                              const List& conv_tract,
                              double mu_total,
                              double hotspot_alpha = 0.0) {
    int k = (int)arr.size();
    if (k == 0) return;

    // === LOCAL DUPLICATION ===
    if (p_local_dup > 0.0) {
        k = (int)arr.size();
        std::vector<int> triggered;
        if (hotspot_alpha > 0.0) {
            // HOTSPOT FEEDBACK (g9): the per-monomer trigger probability is
            // non-uniform, p_i = p_local_dup * k * w_i / sum(w), with
            // w_i = (same-chromosome lineage-family size)^alpha. Singletons get
            // w = 1; a monomer with d additional identical copies gets
            // (1+d)^alpha. MEAN-PRESERVING per chromosome by construction
            // (sum p_i = p_local_dup * k), so alpha redistributes duplication
            // events without adding mass -- the size governor's balance
            // arithmetic is untouched and alpha is not confounded with
            // p_local_dup. (p_i > 1 saturates in the Bernoulli draw; only an
            // extreme-alpha corner ever gets there.) Consumes exactly k
            // uniforms, same as the uniform path.
            std::unordered_map<long long, double> famw;
            famw.reserve(k * 2);
            for (int i = 0; i < k; i++) famw[arr[i].fam] += 1.0;
            for (auto& kv : famw) kv.second = std::pow(kv.second, hotspot_alpha);
            double W = 0.0;
            for (int i = 0; i < k; i++) W += famw[arr[i].fam];
            double scale = (W > 0.0) ? p_local_dup * (double)k / W : 0.0;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < famw[arr[i].fam] * scale)
                    triggered.push_back(i);
            }
        } else {
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_local_dup) triggered.push_back(i);
            }
        }
        for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
            int idx = triggered[ti];
            int ck = (int)arr.size();
            int max_chunk = ck - idx;
            if (max_chunk <= 0) continue;
            int chunk_size = sample_chunk_size(local_dist, max_chunk);
            int start = idx;
            int end = std::min(start + chunk_size - 1, ck - 1);
            for (int j = start; j <= end; j++) materialize(arr[j], mu_total);
            std::vector<Monomer> chunk;
            chunk.reserve(end - start + 1);
            for (int j = start; j <= end; j++) {
                chunk.emplace_back(arr[j].seq, arr[j].dir);
                chunk.back().fam = arr[j].fam;   // copy joins source's family
            }
            arr.insert(arr.begin() + end + 1,
                       std::make_move_iterator(chunk.begin()),
                       std::make_move_iterator(chunk.end()));
        }
    }

    // === DISTAL DUPLICATION (within array) ===
    if (p_distal_dup > 0.0 && R::runif(0.0, 1.0) < p_distal_dup) {
        k = (int)arr.size();
        if (k > 0) {
            int idx = (int)(R::runif(0.0, 1.0) * k);
            if (idx >= k) idx = k - 1;
            int max_chunk = k - idx;
            int chunk_size = sample_chunk_size(distal_dist, max_chunk);
            int start = idx;
            int end = std::min(start + chunk_size - 1, k - 1);
            for (int j = start; j <= end; j++) materialize(arr[j], mu_total);
            bool invert = R::runif(0.0, 1.0) < p_invert_distal;
            std::vector<Monomer> chunk;
            chunk.reserve(end - start + 1);
            if (invert) {
                for (int j = end; j >= start; j--) {
                    chunk.emplace_back(arr[j].seq, -arr[j].dir);
                    chunk.back().fam = arr[j].fam;
                }
            } else {
                for (int j = start; j <= end; j++) {
                    chunk.emplace_back(arr[j].seq, arr[j].dir);
                    chunk.back().fam = arr[j].fam;
                }
            }
            int n_positions = k + 1;
            int insert_at = (int)(R::runif(0.0, 1.0) * n_positions);
            if (insert_at >= n_positions) insert_at = n_positions - 1;
            arr.insert(arr.begin() + insert_at,
                       std::make_move_iterator(chunk.begin()),
                       std::make_move_iterator(chunk.end()));
        }
    }

    // === CHUNK DELETION ===
    if (p_del_chunk > 0.0) {
        k = (int)arr.size();
        std::vector<int> triggered;
        for (int i = 0; i < k; i++) {
            if (R::runif(0.0, 1.0) < p_del_chunk) triggered.push_back(i);
        }
        for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
            int idx = triggered[ti];
            int ck = (int)arr.size();
            if (idx >= ck) continue;
            int max_chunk = ck - idx;
            if (max_chunk <= 0) continue;
            int chunk_size = sample_chunk_size(del_dist, max_chunk);
            int start = idx;
            int end = std::min(start + chunk_size - 1, ck - 1);
            arr.erase(arr.begin() + start, arr.begin() + end + 1);
        }
    }

    // === GENE CONVERSION (homogenization) ===
    if (p_conversion > 0.0) {
        k = (int)arr.size();
        if (k >= 2) {
            std::vector<int> recip;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_conversion) recip.push_back(i);
            }
            for (size_t ri = 0; ri < recip.size(); ri++) {
                int i = recip[ri];
                int off = sample_chunk_size(conv_dist, k - 1);
                int sign = (R::runif(0.0, 1.0) < 0.5) ? -1 : 1;
                int donor = i + sign * off;
                if (donor < 0)      donor = i + off;
                if (donor > k - 1)  donor = i - off;
                if (donor < 0 || donor > k - 1 || donor == i) continue;
                int room = 1 + std::min(k - 1 - i, k - 1 - donor);
                int tract = sample_chunk_size(conv_tract, room);
                for (int j = 0; j < tract; j++) {
                    int rpos = i + j, dpos = donor + j;
                    if (rpos > k - 1 || dpos > k - 1) break;
                    materialize(arr[dpos], mu_total);
                    arr[rpos].seq = arr[dpos].seq;
                    arr[rpos].pending_gens = 0;
                    arr[rpos].fam = arr[dpos].fam;  // recipient joins donor's family
                }
            }
        }
    }

    // === TERMINAL (LONG-RANGE) DELETION ===
    if (p_distal_del > 0.0 && R::runif(0.0, 1.0) < p_distal_del) {
        k = (int)arr.size();
        if (k > 0) {
            int chunk_size = sample_chunk_size(distal_dist, k);
            if (chunk_size > k) chunk_size = k;
            arr.erase(arr.end() - chunk_size, arr.end());
        }
    }
}

// ---------------------------------------------------------------------------
// Core simulation engine
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List sim_core_cpp(IntegerVector ancestor_seq_r,
                  int init_k0,
                  int max_units,
                  int max_t,
                  int hard_cap,
                  double p_local_dup,
                  List local_dist,
                  double p_distal_dup,
                  List distal_dist,
                  double p_invert_distal,
                  double p_del_chunk,
                  List del_dist,
                  double p_distal_del,
                  double p_conversion,
                  List conv_dist,
                  List conv_tract,
                  double mu_total,
                  bool verbose) {

    // Convert ancestor sequence: R 1-indexed -> C++ 0-indexed
    int L = ancestor_seq_r.size();
    std::vector<int> ancestor_seq(L);
    for (int i = 0; i < L; i++) {
        ancestor_seq[i] = ancestor_seq_r[i] - 1;
    }

    // Initialize array with init_k0 identical monomers
    std::vector<Monomer> arr;
    arr.reserve(std::min(max_units, hard_cap) + 1000);
    for (int i = 0; i < init_k0; i++) {
        arr.emplace_back(ancestor_seq, 1);
    }

    // Trajectory: array length per generation
    std::vector<int> l_vec;
    l_vec.reserve(100000);
    l_vec.push_back(init_k0);

    int t = 0;
    bool hit_hard_cap = false;

    while ((int)arr.size() < max_units && t < max_t) {
        t++;
        int k = (int)arr.size();

        if (k > hard_cap) {
            hit_hard_cap = true;
            break;
        }
        // If array went extinct, re-initialize and reset clock
        if (k == 0) {
            arr.clear();
            for (int i = 0; i < init_k0; i++) {
                arr.emplace_back(ancestor_seq, 1);
            }
            t = 0;  // reset clock — new array starts fresh
            continue;
        }

        if (verbose && t % 500 == 0) {
            Rprintf("  Gen %d - units: %d\n", t, k);
        }
        if (t % 100 == 0) Rcpp::checkUserInterrupt();

        // === LOCAL DUPLICATION ===
        if (p_local_dup > 0.0) {
            k = (int)arr.size();
            std::vector<int> triggered;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_local_dup) {
                    triggered.push_back(i);
                }
            }
            for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
                int idx = triggered[ti];
                int ck = (int)arr.size();
                int max_chunk = ck - idx;
                if (max_chunk <= 0) continue;

                int chunk_size = sample_chunk_size(local_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, ck - 1);

                // Materialize source monomers before copying
                for (int j = start; j <= end; j++) {
                    materialize(arr[j], mu_total);
                }

                // Create copies
                std::vector<Monomer> chunk;
                chunk.reserve(end - start + 1);
                for (int j = start; j <= end; j++) {
                    chunk.emplace_back(arr[j].seq, arr[j].dir);
                }

                // Insert right after source (tandem)
                arr.insert(arr.begin() + end + 1,
                           std::make_move_iterator(chunk.begin()),
                           std::make_move_iterator(chunk.end()));
            }
        }

        // === DISTAL DUPLICATION ===
        // Per-ARRAY probability (not per-unit): models chromosome-level events
        // like unequal crossover. At most one event per generation.
        if (p_distal_dup > 0.0 && R::runif(0.0, 1.0) < p_distal_dup) {
            k = (int)arr.size();
            if (k > 0) {
                int idx = (int)(R::runif(0.0, 1.0) * k);
                if (idx >= k) idx = k - 1;
                int max_chunk = k - idx;

                int chunk_size = sample_chunk_size(distal_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, k - 1);

                // Materialize source
                for (int j = start; j <= end; j++) {
                    materialize(arr[j], mu_total);
                }

                // Create copies, possibly inverted
                bool invert = R::runif(0.0, 1.0) < p_invert_distal;
                std::vector<Monomer> chunk;
                chunk.reserve(end - start + 1);
                if (invert) {
                    for (int j = end; j >= start; j--) {
                        chunk.emplace_back(arr[j].seq, -arr[j].dir);
                    }
                } else {
                    for (int j = start; j <= end; j++) {
                        chunk.emplace_back(arr[j].seq, arr[j].dir);
                    }
                }

                // Insert at random position (anywhere in array)
                int n_positions = k + 1;
                int insert_at = (int)(R::runif(0.0, 1.0) * n_positions);
                if (insert_at >= n_positions) insert_at = n_positions - 1;

                arr.insert(arr.begin() + insert_at,
                           std::make_move_iterator(chunk.begin()),
                           std::make_move_iterator(chunk.end()));
            }
        }

        // === CHUNK DELETION ===
        if (p_del_chunk > 0.0) {
            k = (int)arr.size();
            std::vector<int> triggered;
            for (int i = 0; i < k; i++) {
                if (R::runif(0.0, 1.0) < p_del_chunk) {
                    triggered.push_back(i);
                }
            }
            for (int ti = (int)triggered.size() - 1; ti >= 0; ti--) {
                int idx = triggered[ti];
                int ck = (int)arr.size();
                if (idx >= ck) continue;
                int max_chunk = ck - idx;
                if (max_chunk <= 0) continue;

                int chunk_size = sample_chunk_size(del_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, ck - 1);

                arr.erase(arr.begin() + start, arr.begin() + end + 1);
            }
        }

        // === GENE CONVERSION (homogenization) ===
        // Non-allelic gene conversion: a recipient monomer's SEQUENCE is over-
        // written by a NEARBY donor's current sequence (double-strand-break repair
        // templated on a nearby repeat unit). No length change -- this homogenizes
        // (reduces divergence, creates identical copies) WITHOUT adding units.
        // Donor separation ~ conv_dist (local kernel); tract ~ conv_tract (short,
        // typically ~1 monomer). Recipient orientation (dir) is unchanged.
        if (p_conversion > 0.0) {
            k = (int)arr.size();
            if (k >= 2) {
                std::vector<int> recip;
                for (int i = 0; i < k; i++) {
                    if (R::runif(0.0, 1.0) < p_conversion) recip.push_back(i);
                }
                for (size_t ri = 0; ri < recip.size(); ri++) {
                    int i = recip[ri];
                    int off = sample_chunk_size(conv_dist, k - 1);   // donor separation >=1
                    int sign = (R::runif(0.0, 1.0) < 0.5) ? -1 : 1;
                    int donor = i + sign * off;
                    if (donor < 0)      donor = i + off;             // reflect at ends
                    if (donor > k - 1)  donor = i - off;
                    if (donor < 0 || donor > k - 1 || donor == i) continue;
                    int room = 1 + std::min(k - 1 - i, k - 1 - donor);
                    int tract = sample_chunk_size(conv_tract, room);
                    for (int j = 0; j < tract; j++) {
                        int rpos = i + j, dpos = donor + j;
                        if (rpos > k - 1 || dpos > k - 1) break;
                        materialize(arr[dpos], mu_total);   // resolve donor to current seq
                        arr[rpos].seq = arr[dpos].seq;      // overwrite recipient sequence
                        arr[rpos].pending_gens = 0;         // freshly homogenized copy
                    }
                }
            }
        }

        // === TERMINAL (LONG-RANGE) DELETION ===
        // Per-ARRAY probability (mirror of distal duplication): with prob
        // p_distal_del, remove a chunk of size ~distal_dist from the END of the
        // array (telomere-proximal truncation -- arrays lose their distal end).
        if (p_distal_del > 0.0 && R::runif(0.0, 1.0) < p_distal_del) {
            k = (int)arr.size();
            if (k > 0) {
                int chunk_size = sample_chunk_size(distal_dist, k);
                if (chunk_size > k) chunk_size = k;
                arr.erase(arr.end() - chunk_size, arr.end());
            }
        }

        // === INCREMENT PENDING MUTATIONS FOR ALL MONOMERS ===
        for (auto& m : arr) {
            m.pending_gens++;
        }

        l_vec.push_back((int)arr.size());
    }

    // === MATERIALIZE ALL REMAINING MUTATIONS ===
    int n_final = (int)arr.size();
    for (auto& m : arr) {
        materialize(m, mu_total);
    }

    // === CONVERT TO R OUTPUT ===
    const char bases[] = {'A', 'C', 'G', 'T'};
    CharacterVector seq_out(n_final);
    CharacterVector dir_out(n_final);

    for (int i = 0; i < n_final; i++) {
        std::string s(arr[i].seq.size(), ' ');
        for (int j = 0; j < (int)arr[i].seq.size(); j++) {
            s[j] = bases[arr[i].seq[j]];
        }
        seq_out[i] = s;
        dir_out[i] = (arr[i].dir > 0) ? "+" : "-";
    }

    return List::create(
        Named("seqs") = seq_out,
        Named("dirs") = dir_out,
        Named("l_vec") = wrap(l_vec),
        Named("total_gens") = t,
        Named("hit_hard_cap") = hit_hard_cap
    );
}

// ---------------------------------------------------------------------------
// Multi-chromosome (single-haplotype) simulation engine
// ---------------------------------------------------------------------------
// Holds K chromosome arrays in one genome. Each generation:
//   1. for each chromosome: apply the within-array moves (shared helper),
//   2. with prob p_translocation (per genome per generation), copy a contiguous
//      chunk from a random DONOR chromosome into a random RECIPIENT chromosome
//      at a random position, inverted w.p. p_invert_transloc.
// AGE MODE ONLY: every chromosome is seeded AT its target size and the loop
// runs for exactly n_generations of governed steady-state turnover.
//
// One ancestral monomer seeds all chromosomes (shared origin), so without
// translocation, identical sequences shared ACROSS different chromosomes only
// arise from the (vanishing) chance of independent convergence — i.e. ~0.

// Monomers used to re-seed a chromosome in the (practically impossible under
// the governor) event it deletes to extinction.
static const int RESEED_K0 = 10;

// [[Rcpp::export]]
List sim_genome_cpp(IntegerVector ancestor_seq_r,
                    int K,
                    IntegerVector target_sizes_r,
                    IntegerVector hard_caps_r,
                    int n_generations,
                    double size_band,
                    double p_local_dup,
                    List local_dist,
                    double p_distal_dup,
                    List distal_dist,
                    double p_invert_distal,
                    List del_dist,
                    double p_translocation,
                    List transloc_dist,
                    double p_invert_transloc,
                    double mu_total,
                    bool verbose,
                    double p_conversion = 0.0,
                    Nullable<List> conv_dist_ = R_NilValue,
                    Nullable<List> conv_tract_ = R_NilValue,
                    double p_conv_interchrom = 0.0,
                    double hotspot_alpha = 0.0) {

    // Non-allelic gene conversion (homogenization) is OFF by default so this
    // engine stays bit-identical to the released behaviour when unused.
    List conv_dist  = conv_dist_.isNull()
        ? List::create(_["type"] = "gamma", _["shape"] = 1.0, _["scale"] = 20.0)
        : List(conv_dist_);
    List conv_tract = conv_tract_.isNull()
        ? List::create(_["type"] = "gamma", _["shape"] = 1.0, _["scale"] = 1.0)
        : List(conv_tract_);

    // Split the single conversion rate between a WITHIN-chromosome path (donor
    // drawn from conv_dist a short distance away, handled by the shared
    // within-array helper) and an INTER-chromosomal path (donor drawn uniformly
    // from a different chromosome, handled at genome level below). Homogenization
    // by unbiased conversion runs on a coalescent timescale ~ (array size)/rate,
    // so short-range donors homogenize far more effectively than distant ones;
    // the inter-chrom path is therefore NOT for building deep identical families
    // but for seeding identical monomers ACROSS chromosomes without translocating
    // a contiguous block.
    if (p_conv_interchrom < 0.0 || p_conv_interchrom > 1.0)
        stop("p_conv_interchrom must be in [0, 1]");
    double p_conv_within = p_conversion * (1.0 - p_conv_interchrom);
    double p_conv_xchrom = p_conversion * p_conv_interchrom;

    if (K < 1) stop("K must be >= 1");
    if (n_generations < 1) stop("n_generations must be >= 1");
    if ((int)target_sizes_r.size() != K) stop("target_sizes must have length K");
    if ((int)hard_caps_r.size() != K) stop("hard_caps must have length K");

    // Convert ancestor: R 1-indexed -> C++ 0-indexed
    int L = ancestor_seq_r.size();
    std::vector<int> ancestor_seq(L);
    for (int i = 0; i < L; i++) ancestor_seq[i] = ancestor_seq_r[i] - 1;

    std::vector<int> target_sizes(K), hard_caps(K);
    for (int c = 0; c < K; c++) {
        target_sizes[c] = target_sizes_r[c];
        hard_caps[c]    = hard_caps_r[c];
    }

    // Initialize K arrays from the SAME ancestral monomer.
    // Seed each chromosome AT its target size (all identical ancestral
    // monomers, divergence 0). The fixed n_generations of governed turnover
    // then builds divergence (~mu*age) while holding redundancy — a "fully
    // grown young array" aging in place.
    std::vector<std::vector<Monomer> > genome(K);
    for (int c = 0; c < K; c++) {
        int seed_n = target_sizes[c];
        genome[c].reserve(std::min(seed_n, hard_caps[c]) + 1000);
        for (int i = 0; i < seed_n; i++) {
            genome[c].emplace_back(ancestor_seq, 1);
            genome[c].back().fam = (long long)c;  // one seed family per chrom
        }
    }

    // HOTSPOT FEEDBACK (g9) family bookkeeping. next_fam hands out fresh ids;
    // p_fam_break is the exact per-generation probability that >=1 substitution
    // lands anywhere in a monomer (per-base change prob under one JC generation
    // is mu_total), i.e. the rate at which true exact-sequence identity breaks.
    // The whole mechanism is guarded by hotspot_alpha > 0 and consumes ZERO RNG
    // when off, keeping alpha = 0 bit-identical to the released engine.
    long long next_fam = (long long)K;
    double p_fam_break = 1.0 - std::pow(1.0 - mu_total, (double)L);

    int t = 0;
    bool hit_hard_cap = false;
    long n_transloc = 0;

    // STEADY-STATE GENOME DYNAMICS (age-driven).
    // ---------------------------------------------------------------------
    // The loop runs for EXACTLY n_generations generations (age is the
    // divergence dial: with mu fixed, divergence ~= mu * age). Every
    // chromosome is seeded AT its own target size and runs STEADY-STATE
    // TURNOVER -- local duplication keeps firing at p_local_dup (fresh
    // identical twins) while a single GLOBAL (whole-genome) SIZE GOVERNOR
    // sets one genome-wide chunk-deletion prob that matches the TOTAL
    // duplication mass, with a proportional correction that pulls the SUM of
    // all chromosome sizes back toward the total target only when it drifts
    // outside a +/- band deadzone. Individual chromosomes are FREE to
    // expand/shrink (translocation moves mass between them); only the total
    // is held. So redundancy is steady and translocated blocks persist
    // instead of being deleted by a per-chrom governor reacting to recipient
    // overshoot. The governor uses dup/del chunk MEANS to balance per-unit
    // mass: the per-unit chunk deletion absorbs BOTH the local AND the
    // within-array distal duplication mass.
    int  t_limit  = n_generations;
    double band   = (size_band > 0.0) ? size_band : 0.10;  // +/- fraction of TOTAL target (deadband)
    double mean_local  = dist_mean(local_dist);
    double mean_del    = dist_mean(del_dist);
    double mean_distal = dist_mean(distal_dist);
    // proportional governor gain (per-unit deletion-prob correction per unit of
    // fractional size error); chosen so a full-band overshoot adds ~the balance
    // deletion again. Clamped downstream so it can never exceed a sane rate.
    double gov_gain = 0.5;
    // WEAK per-chromosome restoring gain: on top of the global total governor,
    // each chromosome is gently pulled back toward ITS OWN target size when it
    // drifts outside the +/- band, so chromosomes keep their characteristic sizes
    // (the long-term per-chrom size stability seen in the real genomes) instead
    // of random-walking to the genome mean. MUCH weaker than gov_gain so the
    // relaxation time (~1/perchrom_gain gens) is long: a translocation swell
    // persists as a transient bump and erodes only slowly (at the background
    // turnover rate, NOT preferentially), keeping the cross-chrom signal visible.
    // Since the per-chrom targets SUM to the total target, this restoring is
    // consistent with the global governor (they pull the same direction).
    double perchrom_gain = 0.01;

    while (t < t_limit) {
        t++;

        // Extinction guard per chromosome (re-seed) and cap flag (informational
        // only -- does NOT terminate the genome).
        for (int c = 0; c < K; c++) {
            if ((int)genome[c].size() > hard_caps[c]) hit_hard_cap = true;
            if ((int)genome[c].size() == 0) {
                // re-seed an extinct chromosome from ancestor (one new family)
                long long f = next_fam++;
                for (int i = 0; i < RESEED_K0; i++) {
                    genome[c].emplace_back(ancestor_seq, 1);
                    genome[c].back().fam = f;
                }
            }
        }

        if (verbose && t % 500 == 0) {
            Rprintf("  Gen %d:", t);
            for (int c = 0; c < K; c++) Rprintf(" %d", (int)genome[c].size());
            Rprintf("\n");
        }
        if (t % 100 == 0) Rcpp::checkUserInterrupt();

        // === 1. WITHIN-ARRAY MOVES (governed) ===

        // GLOBAL (WHOLE-GENOME) SIZE GOVERNOR.
        // ----------------------------------------------------------------------
        // Size is controlled on the TOTAL monomer count summed over ALL K
        // chromosomes, NOT per chromosome. Individual chromosomes are free to
        // expand and shrink -- a chromosome that RECEIVES a translocated block
        // stays larger and its DONOR stays smaller -- while only the genome-wide
        // SUM is held near the total target (within a +/- `band` deadzone where
        // it floats freely on the net-zero balance, with a proportional pull only
        // once it drifts outside the band).
        //
        // WHY: the old per-chromosome governor raised a RECIPIENT chromosome's own
        // deletion the instant a translocation pushed it above its band, eroding
        // the very cross-chromosome block we are trying to accumulate. A single
        // genome-wide deletion rate lets translocated blocks persist (they only
        // erode at the background turnover rate, like any monomer), so standing
        // cross-chrom sharing -- the fingerprint's off-diagonal ink -- can build
        // up. Chromosomes changing size via translocation becomes an emergent,
        // inferrable feature instead of something the governor cancels.
        //
        // The single genome-wide per-unit deletion prob balances the TOTAL
        // duplication mass (local + within-array distal, summed over chromosomes):
        //   local dup mass  = p_local_dup  * total_sz   * mean_local  (per-unit x all monomers)
        //   distal dup mass = p_distal_dup * mean_distal * K          (per-array event x K chroms)
        //   deletion mass   = p_del        * total_sz   * mean_del    (per-unit x all monomers)
        // => p_del_balance  = p_local_dup*mean_local/mean_del
        //                   + p_distal_dup*mean_distal*K/(total_sz*mean_del)
        // plus gov_gain * (banded fractional total-size error) / mean_del.
        double total_target = 0.0, total_sz = 0.0;
        for (int c = 0; c < K; c++) {
            total_target += (double)target_sizes[c];
            total_sz     += (double)genome[c].size();
        }
        double p_del_global = 0.0;
        if (mean_del > 0.0 && total_sz > 0.0) {
            double p_del_balance = p_local_dup * mean_local / mean_del
                                 + p_distal_dup * mean_distal * (double)K
                                   / (total_sz * mean_del);
            // banded fractional size error: zero inside +/- band (free float),
            // positive above the band (delete more), negative below (delete less).
            double err_raw = (total_target > 0.0)
                             ? (total_sz - total_target) / total_target : 0.0;
            double err = 0.0;
            if      (err_raw >  band) err = err_raw - band;
            else if (err_raw < -band) err = err_raw + band;
            p_del_global = p_del_balance + gov_gain * err / mean_del;
            if (p_del_global < 0.0) p_del_global = 0.0;
            if (p_del_global > 0.5) p_del_global = 0.5;
        }

        for (int c = 0; c < K; c++) {
            int sz = (int)genome[c].size();

            // GLOBAL-governed steady-state turnover with a WEAK per-chrom
            // restoring anchor. Local dup AND within-array distal dup keep
            // firing (fresh twins / off-diagonal copies); the genome-wide
            // governed deletion p_del_global absorbs the total duplication mass,
            // and a weak per-chrom term gently pulls this chromosome back toward
            // its OWN target when it drifts outside the band (long-term size
            // stability without erasing transient translocation swells).
            // A chromosome at/over its hard cap stops GROWING (dup off) but still
            // DELETES, keeping memory bounded while the rest of the genome carries
            // the size.
            double p_local_eff  = p_local_dup;
            double p_distal_eff = p_distal_dup;
            if (sz >= hard_caps[c]) { p_local_eff = 0.0; p_distal_eff = 0.0; }
            // weak per-chrom restoring toward own target (banded deadzone)
            double p_del_c = p_del_global;
            double tgt_c   = (double)target_sizes[c];
            if (tgt_c > 0.0 && mean_del > 0.0) {
                double e_raw = ((double)sz - tgt_c) / tgt_c;
                double e = 0.0;
                if      (e_raw >  band) e = e_raw - band;
                else if (e_raw < -band) e = e_raw + band;
                p_del_c = p_del_global + perchrom_gain * e / mean_del;
                if (p_del_c < 0.0) p_del_c = 0.0;
                if (p_del_c > 0.5) p_del_c = 0.5;
            }
            apply_within_array_moves(genome[c],
                                     p_local_eff, local_dist,
                                     p_distal_eff, distal_dist, p_invert_distal,
                                     p_del_c, del_dist, /*p_distal_del=*/0.0,
                                     p_conv_within, conv_dist, conv_tract,
                                     mu_total, hotspot_alpha);
        }

        // === 2. INTER-CHROMOSOMAL TRANSLOCATION (per genome) ===
        if (K >= 2 && p_translocation > 0.0 &&
            R::runif(0.0, 1.0) < p_translocation) {

            // pick donor and a distinct recipient
            int donor = (int)(R::runif(0.0, 1.0) * K);
            if (donor >= K) donor = K - 1;
            int recip = (int)(R::runif(0.0, 1.0) * K);
            if (recip >= K) recip = K - 1;
            if (recip == donor) recip = (recip + 1) % K;

            int dk = (int)genome[donor].size();
            // A recipient already at/over its hard cap can still DONATE but must
            // not GROW unbounded -- skip insertion into an over-cap recipient.
            if (dk > 0 && (int)genome[recip].size() < hard_caps[recip]) {
                int idx = (int)(R::runif(0.0, 1.0) * dk);
                if (idx >= dk) idx = dk - 1;
                int max_chunk = dk - idx;
                int chunk_size = sample_chunk_size(transloc_dist, max_chunk);
                int start = idx;
                int end = std::min(start + chunk_size - 1, dk - 1);

                // resolve donor monomers to current sequence before copying
                for (int j = start; j <= end; j++)
                    materialize(genome[donor][j], mu_total);

                bool invert = R::runif(0.0, 1.0) < p_invert_transloc;
                std::vector<Monomer> chunk;
                chunk.reserve(end - start + 1);
                if (invert) {
                    for (int j = end; j >= start; j--) {
                        chunk.emplace_back(genome[donor][j].seq,
                                           -genome[donor][j].dir);
                        chunk.back().fam = genome[donor][j].fam;
                    }
                } else {
                    for (int j = start; j <= end; j++) {
                        chunk.emplace_back(genome[donor][j].seq,
                                           genome[donor][j].dir);
                        chunk.back().fam = genome[donor][j].fam;
                    }
                }

                // random insertion position in the recipient
                int rk = (int)genome[recip].size();
                int n_positions = rk + 1;
                int insert_at = (int)(R::runif(0.0, 1.0) * n_positions);
                if (insert_at >= n_positions) insert_at = n_positions - 1;
                genome[recip].insert(genome[recip].begin() + insert_at,
                                     std::make_move_iterator(chunk.begin()),
                                     std::make_move_iterator(chunk.end()));
                n_transloc++;
            }
        }

        // === 3. INTER-CHROMOSOMAL GENE CONVERSION (per genome) ===
        // Each monomer is a recipient w.p. p_conv_xchrom; its donor is a
        // uniformly chosen position on a uniformly chosen DIFFERENT chromosome.
        // Copies sequence only (orientation is a property of the locus, not the
        // converted tract), mirroring the within-array conversion path.
        if (K >= 2 && p_conv_xchrom > 0.0) {
            for (int c = 0; c < K; c++) {
                int kc = (int)genome[c].size();
                if (kc < 1) continue;
                for (int i = 0; i < kc; i++) {
                    if (R::runif(0.0, 1.0) >= p_conv_xchrom) continue;
                    int d = (int)(R::runif(0.0, 1.0) * (K - 1));
                    if (d >= K - 1) d = K - 2;
                    if (d >= c) d++;                       // uniform over chroms != c
                    int kd = (int)genome[d].size();
                    if (kd < 1) continue;
                    int dpos = (int)(R::runif(0.0, 1.0) * kd);
                    if (dpos >= kd) dpos = kd - 1;
                    int room  = 1 + std::min(kc - 1 - i, kd - 1 - dpos);
                    int tract = sample_chunk_size(conv_tract, room);
                    for (int j = 0; j < tract; j++) {
                        int rpos = i + j, dp = dpos + j;
                        if (rpos > kc - 1 || dp > kd - 1) break;
                        materialize(genome[d][dp], mu_total);
                        genome[c][rpos].seq = genome[d][dp].seq;
                        genome[c][rpos].pending_gens = 0;
                        genome[c][rpos].fam = genome[d][dp].fam;
                    }
                    i += tract - 1;                        // skip the converted tract
                }
            }
        }

        // === 4. INCREMENT PENDING MUTATIONS FOR ALL MONOMERS, ALL CHROMS ===
        // With hotspot feedback on, each monomer also leaves its lineage family
        // (fresh singleton id) w.p. p_fam_break -- the exact-identity breakage
        // rate implied by mu_total. Deferred mutation means the ACTUAL
        // substitution lands later (at materialization) with independent RNG;
        // family membership tracks identity in distribution, not per-base
        // realization. Zero RNG when hotspot_alpha == 0.
        for (int c = 0; c < K; c++) {
            for (auto& m : genome[c]) {
                m.pending_gens++;
                if (hotspot_alpha > 0.0 &&
                    R::runif(0.0, 1.0) < p_fam_break) m.fam = next_fam++;
            }
        }
    }

    // === MATERIALIZE ALL REMAINING MUTATIONS ===
    for (int c = 0; c < K; c++)
        for (auto& m : genome[c]) materialize(m, mu_total);

    // === BUILD FLAT OUTPUT: chrom, num, bponly, dir ===
    int n_total = 0;
    for (int c = 0; c < K; c++) n_total += (int)genome[c].size();

    IntegerVector chrom_out(n_total);
    IntegerVector num_out(n_total);
    CharacterVector seq_out(n_total);
    CharacterVector dir_out(n_total);

    const char bases[] = {'A', 'C', 'G', 'T'};
    int row = 0;
    for (int c = 0; c < K; c++) {
        int nc = (int)genome[c].size();
        for (int i = 0; i < nc; i++) {
            std::string s(genome[c][i].seq.size(), ' ');
            for (int j = 0; j < (int)genome[c][i].seq.size(); j++)
                s[j] = bases[genome[c][i].seq[j]];
            chrom_out[row] = c + 1;        // 1..K
            num_out[row]   = i + 1;        // 1..nc within chromosome
            seq_out[row]   = s;
            dir_out[row]   = (genome[c][i].dir > 0) ? "+" : "-";
            row++;
        }
    }

    IntegerVector sizes_out(K);
    for (int c = 0; c < K; c++) sizes_out[c] = (int)genome[c].size();

    return List::create(
        Named("chrom")        = chrom_out,
        Named("num")          = num_out,
        Named("bponly")       = seq_out,
        Named("dir")          = dir_out,
        Named("array_sizes")  = sizes_out,
        Named("total_gens")   = t,
        Named("n_transloc")   = (double)n_transloc,
        Named("hit_hard_cap") = hit_hard_cap
    );
}
