type TooltipInstance = { show(): void; dispose(): void };
type BootstrapWithTooltip = {
    Tooltip: { new (el: Element, opts: object): TooltipInstance; getInstance(el: Element): TooltipInstance | null };
};

export function bsTooltip(node: HTMLElement) {
    let instance: TooltipInstance | null = null;
    function init() {
        const bootstrap = (globalThis as typeof globalThis & { bootstrap?: BootstrapWithTooltip }).bootstrap;
        if (!bootstrap) return;
        instance = bootstrap.Tooltip.getInstance(node) ?? new bootstrap.Tooltip(node, { trigger: "hover" });
        if (node.matches(":hover")) instance.show();
    }
    node.addEventListener("mouseenter", init, { once: true });
    return { destroy() { instance?.dispose(); } };
}
