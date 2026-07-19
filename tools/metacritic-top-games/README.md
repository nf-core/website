# Metacritic Top Games → Obsidian

`metacritic_top_games.py` pulls the top-rated video games of all time from
Metacritic and writes them to a Markdown file formatted as a table, ready to
import into [Obsidian](https://obsidian.md).

## Is there an API?

Metacritic has **no official public API**. Its site is powered by an unofficial
JSON backend (`backend.metacritic.com/finder/...`) that this script calls
directly — it returns clean structured data (title, Metascore, platform,
release date) and is far more reliable than scraping HTML. An `apiKey` query
param is required; the public key from Metacritic's own frontend is used by
default and can be overridden with `--api-key` if it rotates.

If the live backend is unreachable (locked-down networks, sandboxes,
Cloudflare blocks), the script automatically falls back to a small **verified
offline snapshot** of the canonical all-time list so you always get output.

## Usage

```bash
pip install requests          # only needed for live fetching

# Default: games with a Metascore over 80, to metacritic_top_games.md
python metacritic_top_games.py

# Only 90+ "must-play" titles, cap at 50, custom output
python metacritic_top_games.py --min-score 90 --limit 50 --output games.md

# Force the bundled offline snapshot (no network)
python metacritic_top_games.py --offline
```

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--min-score` | `80` | Minimum Metascore to include (parameterizable threshold). |
| `--limit` | `200` | Max games to include (`0` = no limit). |
| `--output` / `-o` | `metacritic_top_games.md` | Output Markdown path. |
| `--api-key` | public key | Override the backend `apiKey`. |
| `--offline` | off | Skip the network, use the bundled snapshot. |

## Output

A Markdown file with YAML frontmatter (Obsidian note properties) and a
GitHub-flavored table: rank, game (linked), Metascore, tier, platform, year.
