import { SearchQuery } from "@components/store";

// Field-scoped search for the listing pages (modules, subworkflows, pipelines).
//
// A query is split into space-separated terms (all of which must match — AND).
// A term may be scoped to a single field with a `field:` prefix, otherwise it
// searches name, description and keywords together (the historical behaviour).
//
//   samtools            -> match "samtools" in name, description or keywords
//   k:STR               -> match "str" as a substring of any keyword
//   k:"STR"             -> match a keyword that is exactly "STR"
//   n:samtools          -> match "samtools" in the name
//   d:sort k:alignment  -> both terms must match, each in its own field
//
// Each field has a short and a long alias (e.g. `k:` / `keywords:` / `topics:`).
// Quotes group a term that contains spaces (e.g. k:"structural variant") and,
// for keyword-scoped terms, switch from substring to exact matching.

type SearchField = "all" | "name" | "description" | "keywords";

/** A clickable example query shown in the SearchHint component. */
export interface SearchExample {
    query: string;
    label: string;
}

export interface SearchTerm {
    field: SearchField;
    /** Already lower-cased so matching never has to normalise it again. */
    value: string;
    /** Only meaningful for keyword-scoped terms: require the whole keyword to match. */
    exact: boolean;
}

const FIELD_ALIASES: Record<string, SearchField> = {
    n: "name",
    name: "name",
    d: "description",
    desc: "description",
    description: "description",
    k: "keywords",
    kw: "keywords",
    keyword: "keywords",
    keywords: "keywords",
    topic: "keywords",
    topics: "keywords",
};

// Optional `prefix:` followed by either a "quoted value" or a bare run of
// non-space chars. Covers `k:"STR"`, `k:STR`, `"STR"` and `STR`.
const TERM_RE = /(?:(\w+):)?(?:"([^"]*)"|(\S+))/g;

export function parseSearchQuery(query: string): SearchTerm[] {
    const terms: SearchTerm[] = [];
    let match: RegExpExecArray | null;
    TERM_RE.lastIndex = 0;
    while ((match = TERM_RE.exec(query)) !== null) {
        const [, prefix, quotedValue, bareValue] = match;
        const rawValue = quotedValue ?? bareValue;
        if (rawValue === undefined || rawValue.trim() === "") {
            continue;
        }

        if (prefix) {
            const field = FIELD_ALIASES[prefix.toLowerCase()];
            if (field) {
                terms.push({ field, value: rawValue.toLowerCase(), exact: quotedValue !== undefined });
                continue;
            }
            // Unknown prefix: treat the whole token as an ordinary search term so
            // that things like "read:me" still behave sensibly.
            terms.push({ field: "all", value: `${prefix}:${rawValue}`.toLowerCase(), exact: false });
            continue;
        }

        terms.push({ field: "all", value: rawValue.toLowerCase(), exact: false });
    }
    return terms;
}

interface SearchFields {
    name?: string;
    description?: string;
    /** Any tag-like field: module/subworkflow keywords, pipeline topics, … */
    keywords?: string[];
}

export function matchesSearch(fields: SearchFields, terms: SearchTerm[]): boolean {
    if (terms.length === 0) {
        return true;
    }
    const name = (fields.name || "").toLowerCase();
    const description = (fields.description || "").toLowerCase();
    const keywords = fields.keywords || [];

    // Lower-cased lazily so field-scoped queries (e.g. `n:samtools`) never touch keywords.
    const matchesKeyword = (needle: string, exact: boolean) =>
        keywords.some((keyword) => {
            const lower = keyword.toLowerCase();
            return exact ? lower === needle : lower.includes(needle);
        });

    return terms.every((term) => {
        const needle = term.value;
        switch (term.field) {
            case "name":
                return name.includes(needle);
            case "description":
                return description.includes(needle);
            case "keywords":
                return matchesKeyword(needle, term.exact);
            default:
                return name.includes(needle) || description.includes(needle) || matchesKeyword(needle, false);
        }
    });
}

/** Set the search box to an exact-keyword query — the handler for a clicked keyword/topic badge. */
export function filterByKeyword(keyword: string): void {
    SearchQuery.set(`k:"${keyword}"`);
}
