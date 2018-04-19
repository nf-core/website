//
// nf-core.js
// Custom javascript for the nf-core website
//

$(function () {
    // Enable tooltips
    $('[data-toggle="tooltip"]').tooltip()

    // Enable code highlighting
    hljs.initHighlightingOnLoad();
    // Don't try to guess markdown language to highlight (gets it wrong most of the time)
    hljs.configure({languages: []});
})
