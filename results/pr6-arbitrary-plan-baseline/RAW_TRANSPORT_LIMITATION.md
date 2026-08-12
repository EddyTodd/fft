# Raw transport note

The formal v6 benchmark produced three canonical gzip-compressed CSV session files containing 4,560 observations in total. Their canonical byte counts and SHA-256 digests are recorded in `metadata.json` and are the source of truth for the formal analysis.

This PR was authored through a GitHub connector that accepts UTF-8 contents but does not expose arbitrary binary upload. Connector testing also showed that oversized text payloads can be corrupted in transport. For that reason, only transport shards that could be verified byte-for-byte against their expected Git blob hashes are retained under `raw/`; an incomplete or unverified session is deliberately omitted rather than presented as canonical evidence.

Canonical raw streams:

| Session | gzip bytes | SHA-256 |
|---:|---:|---|
| 0 | 17,926 | `7c04cc342876316101a886f05e9400ad3ad97237553140482326795c07a81b99` |
| 1 | 17,851 | `651d0faf533136e20dd948dd74f65b6ced973b562a668973990ab254cb5d9c47` |
| 2 | 18,059 | `a24fbf740dad4975956f55242132164658785773f32ff2fe64e28652c931508f` |

`ANALYSIS.md` was generated from the complete canonical three-session corpus before publication. The omission here is a repository-transport limitation, not a replacement of raw measurements by aggregates.

A future binary-capable publication pass should attach the three canonical gzip streams and verify these SHA-256 digests before removing this note. Until then, claims in this result package should be interpreted with this provenance limitation explicitly in view.
