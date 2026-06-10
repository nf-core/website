// Satteri replacement for `rehype-pretty-code`.
// Highlights fenced code blocks with Shiki directly (async visitors are supported —
// note this makes markdownToHtml return a Promise) and emits the same hooks the
// site CSS targets: figure[data-rehype-pretty-code-figure], pre/code[data-language]
// [data-theme], span[data-line], and figcaption[data-rehype-pretty-code-title] for
// ```lang title="..." blocks.
// Matches the old options: defaultLang "plaintext", keepBackground (shiki default),
// dual github-dark-dimmed/github-light themes, transformerNotationDiff.
import { transformerNotationDiff } from "@shikijs/transformers";
import { h } from "hastscript";
import { defineHastPlugin } from "satteri";
import { codeToHast } from "shiki";

const THEMES = { dark: "github-dark-dimmed", light: "github-light" };
const THEME_ATTR = Object.values(THEMES).join(" ");

const prettyCodeCompatTransformer = (lang: string) => ({
    name: "pretty-code-compat",
    pre(node: any) {
        node.properties["data-language"] = lang;
        node.properties["data-theme"] = THEME_ATTR;
    },
    code(node: any) {
        node.properties["data-language"] = lang;
        node.properties["data-theme"] = THEME_ATTR;
    },
    line(node: any) {
        node.properties["data-line"] = "";
    },
});

export const prettyCodePlugin = defineHastPlugin({
    name: "pretty-code",
    element: {
        filter: ["pre"],
        async visit(node: any, ctx) {
            const code = node.children?.find((c: any) => c?.type === "element" && c.tagName === "code");
            if (!code) return;
            // Don't re-process blocks we already highlighted (e.g. after wrapNode).
            if (node.properties?.["data-language"] || node.properties?.dataLanguage) return;

            // Satteri keeps the fenced info-string language on data.lang (and as a
            // language-* class, kept as fallback) and the rest on data.meta.
            const classes = Array.isArray(code.properties?.className)
                ? code.properties.className
                : typeof code.properties?.className === "string"
                  ? code.properties.className.split(" ")
                  : [];
            const lang =
                code.data?.lang ??
                classes
                    .find((c: any) => String(c).startsWith("language-"))
                    ?.slice("language-".length) ??
                "plaintext";
            const meta: string = code.data?.meta ?? "";
            const value = ctx.textContent(code).replace(/\n$/, "");

            const options = (language: string) => ({
                lang: language,
                themes: THEMES,
                // Emit only --shiki-light/--shiki-dark variables (no inline color),
                // like rehype-pretty-code — main.scss switches the variables per theme,
                // and an inline color would override the stylesheet.
                defaultColor: false as const,
                cssVariablePrefix: "--shiki-",
                meta: { __raw: meta },
                transformers: [transformerNotationDiff(), prettyCodeCompatTransformer(language)],
            });

            let root: any;
            try {
                root = await codeToHast(value, options(lang));
            } catch {
                root = await codeToHast(value, options("plaintext"));
            }
            const pre = root.children.find((c: any) => c?.type === "element" && c.tagName === "pre");
            if (!pre) return;

            const title = meta.match(/title="([^"]+)"/)?.[1];
            return h("figure", { "data-rehype-pretty-code-figure": "" }, [
                ...(title
                    ? [h("figcaption", { "data-rehype-pretty-code-title": "", "data-theme": THEME_ATTR }, title)]
                    : []),
                pre,
            ]);
        },
    },
});

export default prettyCodePlugin;
