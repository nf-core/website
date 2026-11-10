// Dead-code rules only: formatting is handled by prettier.
import astro from "eslint-plugin-astro";
import svelte from "eslint-plugin-svelte";
import tseslint from "typescript-eslint";

export default [
  { ignores: ["**/dist/**", "**/.astro/**", "**/.netlify/**", "**/node_modules/**"] },
  ...astro.configs["flat/base"],
  ...svelte.configs["flat/base"],
  {
    files: ["**/*.svelte"],
    languageOptions: {
      parserOptions: { parser: tseslint.parser, extraFileExtensions: [".svelte"] },
    },
  },
  {
    files: ["**/*.ts", "**/*.mts"],
    languageOptions: { parser: tseslint.parser },
  },
  {
    // also matches **/*.astro/*.{js,ts}, the <script> blocks extracted by the astro processor
    files: ["**/*.{js,mjs,cjs,ts,mts}", "**/*.astro", "**/*.svelte"],
    plugins: { "@typescript-eslint": tseslint.plugin },
    rules: {
      "@typescript-eslint/no-unused-vars": ["error", { argsIgnorePattern: "^_", varsIgnorePattern: "^_" }],
      "@typescript-eslint/no-unused-expressions": "error",
      "@typescript-eslint/no-useless-constructor": "error",
      "no-unused-labels": "error",
      "no-unused-private-class-members": "error",
      "no-unreachable": "error",
      "no-constant-condition": "error",
      "no-dupe-else-if": "error",
      "no-lone-blocks": "error",
      "no-self-assign": "error",
      "no-useless-catch": "error",
      "no-useless-computed-key": "error",
      "no-useless-concat": "error",
      "no-useless-escape": "error",
      "no-useless-rename": "error",
      "no-useless-return": "error",
    },
  },
  {
    files: ["**/*.astro"],
    rules: {
      // astro/no-unused-css-selector is left off: it can't see classes added at runtime, <style is:inline>
      // or styles forwarded to child components. Run it ad hoc with
      // npx eslint --rule '{"astro/no-unused-css-selector": "error"}' sites
      "astro/no-unused-define-vars-in-style": "error",
    },
  },
  {
    files: ["**/*.svelte"],
    rules: {
      "svelte/no-unused-svelte-ignore": "error",
      "svelte/no-unused-props": "error",
      "svelte/no-dupe-else-if-blocks": "error",
      "svelte/no-useless-mustaches": "error",
      "svelte/no-useless-children-snippet": "error",
    },
  },
];
