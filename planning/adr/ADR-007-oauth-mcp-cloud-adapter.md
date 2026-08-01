# ADR-007 — Authenticated stateless OAuth MCP cloud adapter

- **Status:** accepted
- **Date:** 2026-08-01
- **Issue:** GEN-015

## Context

The six deterministic public operations and the review-only feedback transformer
already exist as a dependency-free Node ESM package. ChatGPT Web and other remote
MCP clients need an authenticated transport, but GEN-050 persistence and GEN-051
stored-record RBAC/consent gates are not dependency-ready and must not be
implemented implicitly.

## Decision

Expose exactly seven fixed read-only MCP tools from a dedicated Supabase project.
Supabase OAuth 2.1 authenticates the existing allowed user; a custom access-token
hook binds the token to the exact MCP resource. The Edge function rejects missing,
wrong-audience, disallowed-user, or stale-contract tokens before dispatch.

The Edge runtime vendors a checksum-verified snapshot of the Node core and must
match its golden outputs. It adds no application table, bucket, queue, result
store, generic dispatch tool, welfare override, or publication path. The seventh
tool only prepares feedback with server-owned package integrity and always returns
`pending_human_review`, `auto_apply:false`, and `target_files:[]`.

A Windows-first Node updater validates the complete source/release matrix,
classifies sensitive changes, requires the pinned owner approval for those
changes, deploys an immutable candidate, verifies live OAuth/MCP behavior, and
restores the previous snapshot on failure. Watch mode debounces changes and uses
the same fail-closed one-shot pipeline; it does not grant AutoResearch mutation.

The OAuth consent UI is hosted by a dedicated Cloudflare Worker because the
Supabase shared Edge domain serves HTML under a sandboxed content type. The
Worker is presentation-only and stores no biological or OAuth result data.

## Consequences

- GEN-013 and GEN-014 remain the only prerequisites of GEN-015.
- GEN-050 and GEN-051 remain unimplemented; cloud transport grants no database,
  stored-record RBAC, or welfare-override authority.
- ChatGPT Web uses its private custom MCP App. The local universal plugin only
  references that registered connection and bundles the verified global skill.
- Public npm-registry publication and public plugin-directory submission remain
  outside this release and require separate owner approval.

## Verification

Completion requires exact seven-tool discovery, OAuth audience/user checks,
Edge/source parity, changed-deploy and unchanged-no-op evidence, rollback evidence,
an installed validated personal plugin, and no unresolved P0/P1 security or
scientific-boundary finding.
