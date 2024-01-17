// taken from https://github.com/JuanM04/portfolio/blob/983b0ed0eabdac37bf8b7912d3e8128a443192b9/src/plugins/mermaid.ts
import type { RemarkPlugin } from '@astrojs/markdown-remark';
import dedent from 'ts-dedent';
import { visit } from 'unist-util-visit';

const escapeMap: Record<string, string> = {
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;',
};

const escapeHtml = (str: string) => str.replace(/[&<>"']/g, (c) => escapeMap[c]);

export const mermaid: RemarkPlugin<[]> = () => (tree) => {
    visit(tree, 'code', (node) => {
        if (node.lang !== 'mermaid') return;

        // @ts-ignore
        node.type = 'html';
        node.value = dedent`
      <div class="mermaid" data-content="${escapeHtml(node.value)}">
        <i class="mt-5 m-auto text-success fa-regular fa-spinner-third fa-spin fa-3x"></i>
        <p>Loading graph</p>
      </div>
    `;
    });
};
