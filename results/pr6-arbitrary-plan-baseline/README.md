# PR6 arbitrary-length planning baseline

This result package is the first formal reusable-plan study for **non-power-of-two complex FFTs** in the repository.

## What is compared

At each prime N in `17,31,61,127,257,509,1009,4093`, the execution matrix contains:

- legacy Bluestein (setup/allocation inside every call);
- planned Bluestein (persistent chirp and convolution-kernel spectrum);
- legacy Rader (setup/allocation inside every call);
- planned Rader (persistent permutations and convolution-kernel spectrum);
- FFTW 3.3.10 `ESTIMATE` persistent plan;
- FFTW 3.3.10 `MEASURE` persistent plan.

Legacy modes are historical API baselines. **Algorithm conclusions use planned-vs-planned comparisons.**

Setup is a separate response variable for planned Bluestein, planned Rader, FFTW `ESTIMATE`, and cold FFTW `MEASURE`.

## Corpus

- 3 independently randomized sessions;
- 8 prime sizes;
- 31 samples per execution mode per size/session;
- 1 setup observation per plan policy per size/session;
- 5 warmups;
- 2 ms calibration target;
- **4,560 raw observations** total;
- every observation retained in deterministic gzip shards.

The exact benchmark source is commit `dc2729623a8b90f97ddad36abc67b705a7b95bc2`. `metadata.json` records source blobs, binary hash, compiler/build, FFTW runtime, raw hashes, validation, and the session-orchestration note.

## Central findings

Planning is not a small wrapper optimization here. Across the formal matrix:

- planned Bluestein is about **3.09–4.47x faster** than legacy Bluestein;
- planned Rader is about **2.35–7.61x faster** than legacy Rader;
- plan setup generally repays itself in roughly **0.8–3.3 transforms**.

The Rader/Bluestein result depends strongly on the structure of `p-1`:

- N=17: planned Rader is **4.34x faster** because its cyclic convolution is length 16 versus Bluestein length 64;
- N=257: planned Rader is **4.33x faster** because its cyclic convolution is length 256 versus Bluestein length 1024;
- when both reductions use the same convolution length, differences are much smaller: Rader wins N=31 by ~1.8%, Bluestein wins N=61,127,509,1009 by roughly 0.8–1.6%, and N=4093 is statistically unresolved.

This evidence rejects the previous blanket rule **prime -> Rader** for the reusable planner. The post-evidence `ArbitraryPlan::Auto` policy therefore uses Rader only when it produces a strictly shorter convolution; otherwise it uses Bluestein. Explicit policies remain available.

FFTW remains faster across the formal matrix. The best local planned reduction ranges from about **1.28x slower** at N=257 to about **3.65x slower** at N=4093 versus FFTW `MEASURE`.

## Reproduce analysis

```bash
python3 tools/analyze_arbitrary.py results/pr6-arbitrary-plan-baseline/raw \
  --seed 20260812 --bootstrap 5000
```

## Scope

These measurements come from a virtualized/containerized AMD EPYC environment. They support implementation-specific planning/crossover conclusions and the structural dispatch rule; they do not justify a universal prime FFT ranking across hardware.
