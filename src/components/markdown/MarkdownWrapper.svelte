<script lang="ts">
    import { currentHeading } from '@components/store';
    import * as icons from 'file-icons-js';
    import 'file-icons-js/css/style.css';
    import { onMount } from 'svelte';
    import mermaid from 'mermaid';

    export let headings: { text: string; slug: string; depth: number; fa_icon?: string }[] = [];
    // find current heading in viewport with IntersectionObserver
    onMount(() => {
        async function renderDiagrams(graphs) {
            mermaid.initialize({
                startOnLoad: false,
                fontFamily: 'var(--sans-font)',
                // @ts-ignore This works, but TS expects a enum for some reason
                theme: window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'default',
            });

            for (const graph of graphs) {
                const content = graph.getAttribute('data-content');
                if (!content) continue;
                let svg = document.createElement('svg');
                const id = (svg.id = 'mermaid-' + Math.round(Math.random() * 100000));
                graph.appendChild(svg);
                mermaid.render(id, content).then((result) => {
                    graph.innerHTML = result.svg;
                });
            }
        }
        const graphs = document.getElementsByClassName('mermaid');
        if (document.getElementsByClassName('mermaid').length > 0) {
            renderDiagrams(graphs);
        }

        const observer = new IntersectionObserver(
            (entries) => {
                entries.forEach((entry) => {
                    if (entry.isIntersecting) {
                        currentHeading.set(entry.target.id);
                    }
                });
            },
            {
                rootMargin: '0px 0px -92% 0px',
            },
        );
        headings.forEach((heading) => {
            const element = document.querySelector('#' + heading.slug);
            observer.observe(element);
        });

        // Add "Copy code" button in code blocks
        const copyButtonLabel = "<i class='fa-regular fa-clipboard'></i>";
        const copiedButtonLabel = `<span class='font-sans-serif'><i class='fa-regular fa-clipboard-check me-2 '></i> Copied</span>`;
        document
            .querySelectorAll("div[data-rehype-pretty-code-fragment] pre:not([data-language='console'])")
            .forEach((block) => {
                // only add button if browser supports Clipboard API
                if (navigator.clipboard) {
                    let button = document.createElement('button');
                    button.classList.add(
                        'copy-code-button',
                        'btn',
                        'btn-sm',
                        'btn-outline-secondary',
                        'position-absolute',
                        'top-0',
                        'end-0',
                        'opacity-50',
                    );
                    button.innerHTML = copyButtonLabel;
                    button.title = 'Copy to clipboard';
                    // add data-bs-toggle="tooltip" to enable Bootstrap tooltips
                    button.setAttribute('data-bs-toggle', 'tooltip');
                    block.classList.add('position-relative');
                    block.appendChild(button);
                    button.addEventListener('click', async (e) => {
                        await copyCode(block, e.currentTarget);
                    });
                }
            });
        async function copyCode(block, button) {
            let code = block.querySelector('code').innerText;
            await navigator.clipboard.writeText(code);
            // visual feedback that task is completed
            button.innerHTML = copiedButtonLabel;
            button.classList.replace('btn-outline-secondary', 'btn-outline-success');
            button.classList.remove('opacity-50');
            setTimeout(() => {
                button.innerHTML = copyButtonLabel;
                button.classList.replace('btn-outline-success', 'btn-outline-secondary');
                button.classList.add('opacity-50');
            }, 1500);
        }
        // Add file icon to code block titles
        document.querySelectorAll('div[data-rehype-pretty-code-title]').forEach((block) => {
            const title = block.textContent;
            const fileIcon = icons.getClass(title);
            let icon: HTMLElement;
            if (fileIcon) {
                icon = document.createElement('span');
                icon.classList.add('mx-2', fileIcon);
            } else {
                icon = document.createElement('i');
                icon.classList.add('fa-regular', 'fa-file-code', 'mx-2');
            }
            block.prepend(icon);
        });
    });
</script>

<div>
    <slot />
</div>
