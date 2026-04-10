const DOCSEARCH_OPTS = {
    appId: "01BY3A8NRJ",
    indexName: "nf-co",
    apiKey: "c726615ab69f88b4e26bc891c97c1808",
    placeholder: "Search",
} as const;

async function mountAll(): Promise<void> {
    const htmlDocument = globalThis["document"] as Document | undefined;
    if (!htmlDocument) return;

    const roots = htmlDocument.querySelectorAll("[data-search-container]");
    if (roots.length === 0) return;

    const mod = await import("@docsearch/js");
    const runDocsearch = mod.default;
    for (const node of Array.from(roots)) {
        if (!(node instanceof HTMLElement) || node.dataset.docsearchReady === "true") continue;
        runDocsearch({ container: node, ...DOCSEARCH_OPTS });
        node.dataset.docsearchReady = "true";
    }
}

export function initDocsearch(): void {
    const htmlDocument = globalThis["document"] as Document | undefined;
    if (!htmlDocument) return;
    const run = () => void mountAll();
    if (htmlDocument?.readyState === "loading") {
        htmlDocument.addEventListener("DOMContentLoaded", run, { once: true });
    } else {
        run();
    }
}
