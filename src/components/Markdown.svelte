<script>
    import Markdown from 'svelte-exmarkdown';
    import emoji from 'remark-emoji';
    import remarkGfm from 'remark-gfm';
    import remarkDirective from 'remark-directive';
    import calloutsPlugin from '/bin/remark-callouts.js';

    export let md;

    // replace > ⚠️ with :::warning amd add ::: at the end
    // md = md.replace(/(> ⚠️)(.*)(\n|\.$)/g, ':::note \n $2 \n:::\n');
    if (md) {
        // replace newline with <br>
        md = md.replace(/(\n)/g, '  \n');
        // escape *
        md = md.replace(/(\*.)/g, '\\$1');
    }
</script>

<Markdown {md} plugins={[{ remarkPlugin: [emoji, remarkGfm, remarkDirective, calloutsPlugin] }]} />
