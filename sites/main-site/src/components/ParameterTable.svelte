<script lang="ts">
    import { onMount } from "svelte";

    type ParamVariant = { type: string | null; pipelines: { name: string; version: string }[] };
    type ParamsData = Record<string, ParamVariant[]>;

    let paramsData = $state<ParamsData | null>(null);
    let loading = $state(true);
    let error = $state<string | null>(null);
    let minOverlap = $state(5);
    let sortParams = $state<"count" | "alpha">("count");
    let sortPipelines = $state<"count" | "alpha">("count");
    let searchQuery = $state("");

    function bsTooltip(node: HTMLElement) {
        let instance: { show(): void; dispose(): void } | null = null;
        function init() {
            const bootstrap = (
                globalThis as typeof globalThis & {
                    bootstrap?: { Tooltip: new (el: Element, opts: object) => { show(): void; dispose(): void } };
                }
            ).bootstrap;
            if (!bootstrap) return;
            instance = new bootstrap.Tooltip(node, { trigger: "hover", placement: "top" });
            if (node.matches(":hover")) instance.show();
        }
        node.addEventListener("mouseenter", init, { once: true });
        return {
            destroy() {
                instance?.dispose();
            },
        };
    }

    onMount(async () => {
        try {
            const res = await fetch("/params.json");
            if (!res.ok) throw new Error(`HTTP ${res.status}`);
            paramsData = await res.json();
        } catch (e) {
            error = String(e);
        } finally {
            loading = false;
        }
    });

    // Map of pipeline name → version string (from any param entry)
    const pipelineVersions = $derived.by(() => {
        const map = new Map<string, string>();
        if (!paramsData) return map;
        for (const variants of Object.values(paramsData))
            for (const { pipelines } of variants)
                for (const { name, version } of pipelines) if (!map.has(name)) map.set(name, version);
        return map;
    });

    const totalPipelines = $derived.by(() => {
        if (!paramsData) return 0;
        const s = new Set<string>();
        for (const variants of Object.values(paramsData))
            for (const { pipelines } of variants) for (const { name } of pipelines) s.add(name);
        return s.size;
    });

    const filteredParams = $derived.by(() => {
        if (!paramsData) return [];
        const maxOv = totalPipelines - 1;
        const q = searchQuery.toLowerCase();

        // Flatten name × type variants into individual rows
        const rows = Object.entries(paramsData).flatMap(([name, variants]) =>
            variants.map(({ type, pipelines }) => ({ name, type, pipelines })),
        );

        let filtered = rows.filter(({ name, pipelines }) => {
            const n = pipelines.length;
            return n >= minOverlap && n <= maxOv && (!q || name.toLowerCase().includes(q));
        });

        filtered.sort((a, b) =>
            sortParams === "count"
                ? b.pipelines.length - a.pipelines.length
                : a.name.localeCompare(b.name) || (a.type ?? "").localeCompare(b.type ?? ""),
        );

        return filtered.map(({ name, type, pipelines }) => ({
            name,
            type,
            label: type ? `${name}<code class="small text-secondary bg-body">${type}</code>` : name,
            count: pipelines.length,
            pipelineSet: new Set(pipelines.map((e) => e.name)),
        }));
    });

    const sortedPipelines = $derived.by(() => {
        if (!paramsData) return [];
        const all = new Set<string>();
        for (const variants of Object.values(paramsData))
            for (const { pipelines } of variants) for (const { name } of pipelines) all.add(name);
        const pipes = [...all];
        if (sortPipelines === "count") {
            const counts = new Map(pipes.map((p) => [p, filteredParams.filter((fp) => fp.pipelineSet.has(p)).length]));
            pipes.sort((a, b) => (counts.get(b) ?? 0) - (counts.get(a) ?? 0));
        } else {
            pipes.sort((a, b) => a.localeCompare(b));
        }
        return pipes;
    });
</script>

{#if loading}
    <div class="d-flex justify-content-center align-items-center py-5">
        <i class="fa-regular fa-spinner-third fa-spin fa-2x text-success me-3"></i>
        <span>Loading params data…</span>
    </div>
{:else if error}
    <div class="alert alert-danger">Failed to load params.json: {error}</div>
{:else}
    <div class="d-flex flex-wrap gap-3 align-items-center mb-3">
        <div class="d-flex align-items-center gap-2">
            <label class="form-label mb-0 small text-nowrap" for="min-overlap">
                Min pipelines: <strong>{minOverlap}</strong>
            </label>
            <input
                id="min-overlap"
                type="range"
                class="form-range"
                min="1"
                max="100"
                bind:value={minOverlap}
                style="width:120px"
            />
        </div>

        <div class="d-flex align-items-center gap-2">
            <label class="form-label mb-0 small text-nowrap" for="sort-params">Sort params:</label>
            <select id="sort-params" class="form-select form-select-sm" style="width:auto" bind:value={sortParams}>
                <option value="count">by pipeline count</option>
                <option value="alpha">alphabetical</option>
            </select>
        </div>

        <div class="d-flex align-items-center gap-2">
            <label class="form-label mb-0 small text-nowrap" for="sort-pipelines">Sort pipelines:</label>
            <select
                id="sort-pipelines"
                class="form-select form-select-sm"
                style="width:auto"
                bind:value={sortPipelines}
            >
                <option value="count">by param count</option>
                <option value="alpha">alphabetical</option>
            </select>
        </div>

        <div style="flex:1;min-width:160px">
            <input
                type="search"
                class="form-control form-control-sm"
                placeholder="Filter params…"
                bind:value={searchQuery}
            />
        </div>

        <small class="text-muted text-nowrap">
            {filteredParams.length} params and {sortedPipelines.length} pipelines
        </small>
    </div>

    <div class="overflow-auto">
        <table class="table table-sm table-bordered table-hover mb-0" role="grid">
            <thead class="sticky-top">
                <tr>
                    <th class="param-col"></th>
                    {#each sortedPipelines as pipeline}
                        {@const usedCount = filteredParams.filter((p) => p.pipelineSet.has(pipeline)).length}
                        {@const version = pipelineVersions.get(pipeline)}
                        {@const isDev = version === "dev"}
                        <th
                            class="p-0 align-bottom text-nowrap pipeline-th"
                            use:bsTooltip
                            data-bs-title="{pipeline} ({version}) · {usedCount} of {filteredParams.length} params"
                        >
                            <div class="pipeline-label fw-normal">
                                <a href="/{pipeline}" class="text-body text-decoration-none pt-2">{pipeline}</a>
                                {#if isDev}
                                    <i class="fa-solid fa-xs fa-fw fa-wrench text-warning" title="dev"></i>
                                {/if}
                            </div>
                        </th>
                    {/each}
                </tr>
            </thead>
            <tbody>
                {#each filteredParams as param}
                    <tr>
                        <td
                            class="param-col text-nowrap small fw-medium py-0 px-2"
                            use:bsTooltip
                            data-bs-html="true"
                            data-bs-title="{param.label} · {param.count} pipeline{param.count !== 1 ? 's' : ''}"
                            >{@html param.label}</td
                        >
                        {#each sortedPipelines as pipeline}
                            <td
                                class="p-0 align-middle cell"
                                class:bg-success={param.pipelineSet.has(pipeline)}
                                use:bsTooltip
                                data-bs-html="true"
                                data-bs-title="{pipeline}:<br/>{param.label}"
                            ></td>
                        {/each}
                    </tr>
                {/each}
            </tbody>
        </table>
    </div>
{/if}

<style>
    .param-col {
        min-width: 12rem;
        position: sticky;
        left: 0;
        background: var(--bs-body-bg);
    }
    thead .param-col {
        z-index: 3;
    }
    .pipeline-label {
        writing-mode: vertical-rl;
        transform: rotate(180deg);
        max-height: 12rem;
        overflow: hidden;
    }
    .cell {
        width: 1rem;
        min-width: 1rem;
    }
</style>
