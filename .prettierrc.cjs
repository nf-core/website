module.exports = {
    // TODO: REMOVE SINGLE QUOTES CONFIG LATER
    // Leaving it for now to minimise diffs, but double quote is more consistent with the rest of nf-core
    singleQuote: true,
    //////

    printWidth: 120,

    // pnpm support, for VSCode?
    plugins: [require.resolve('prettier-plugin-astro')],
    overrides: [
        {
            files: '*.astro',
            options: {
                parser: 'astro',
            },
        },
    ],
};
