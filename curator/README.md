# CHD Atlas curator plugin

Skills that drive agents through the per-gene literature pipeline. They are plain
Markdown: any agent or person can read `skills/*/SKILL.md` and follow it by hand.
`curator-rules.md` holds the house rules every skill loads.

## Load it

From the repository root:

```bash
claude --plugin-dir ./curator
```

That loads the plugin for one session; Claude Code finds `skills/` on its own. There is no
marketplace file, so there is no persistent install.

## Run it

Use `.venv/bin/chd-atlas` from the repository root. For one published gene:

1. `literature-search` runs `chd-atlas literature search <SYMBOL>` and records every hit.
2. `literature-triage` (a cheap model) marks each hit included or excluded in
   `curation/literature/<SYMBOL>.yaml`.
3. `literature-extract` (a mid-tier model) writes findings and adds each paper with
   `chd-atlas literature add-publication`.
4. `literature-verify` (a different model, fresh context) extracts blind into
   `local/literature/<SYMBOL>.verify.yaml`; `chd-atlas literature diff` compares and, with
   `--mark-verified`, records agreement.
5. You review: `.venv/bin/chd-atlas literature review <SYMBOL>` lists the unapproved
   findings with their PubMed links and locators. Approve with `--approve <id> --by "<your
   name>"`, or mark one with `--pending-review <id> --note "<what differs>"`. Only a named
   person approves; the command refuses a name starting `agent:`.

`NCBI_API_KEY` is optional, read from the environment only, and never written to a file.
The search and abstracts commands need the network.

Run `.venv/bin/chd-atlas validate` after each step. The plugin lives outside `curation/`,
so it does not change validation.
