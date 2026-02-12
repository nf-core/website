module.exports = {
    trailingComma: "all",
    printWidth: 120,

    plugins: ["prettier-plugin-astro", "prettier-plugin-svelte"],
    overrides: [
        {
            files: "**/*.astro",
            options: {
                parser: "astro",
            },
        },
        {
            files: "**/*.svelte",
            options: {
                parser: "svelte",
            },
        },
    ],
};
