/** @type {import("prettier").Config} */
export default {
    // TODO: REMOVE SINGLE QUOTES CONFIG LATER
    // Leaving it for now to minimise diffs, but double quote is more consistent with the rest of nf-core
    singleQuote: true,
    trailingComma: 'all',
    //////

    printWidth: 120,

    plugins: ['prettier-plugin-astro', 'prettier-plugin-svelte'],
    overrides: [
        {
            files: '*.astro',
            options: {
                parser: 'astro',
            },
        },
        {
            files: '*.svelte',
            options: {
                parser: 'svelte',
            },
        },
    ],
};
