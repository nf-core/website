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

    // Override the .contains() filter to be case insenstive
    $.expr[":"].contains = $.expr.createPseudo(function(arg) {
        return function( elem ) {
            return $(elem).text().toUpperCase().indexOf(arg.toUpperCase()) >= 0;
        };
    });

    // Filter pipelines with text
    $('.pipelines-toolbar .pipeline-filters input').keyup(function(){
        var ftext = $('.pipelines-toolbar .pipeline-filters input').val();
        $('.pipelines-container .pipeline:contains("'+ftext+'")').show();
        $('.pipelines-container .pipeline:not(:contains("'+ftext+'"))').hide();
        if($('.pipelines-container .pipeline:visible').length == 0){ $('.no-pipelines').show(); }
        else { $('.no-pipelines').hide(); }
    });
    // Filter pipelines with buttons
    $('.pipelines-toolbar .pipeline-filters button').click(function(){
        $(this).blur().toggleClass('active');
        var showclasses = [];
        $('.pipelines-toolbar .pipeline-filters button.active').each(function(){
            showclasses.push($(this).data('target'));
        });
        $('.pipelines-container .pipeline').filter(showclasses.join(', ')).show();
        $('.pipelines-container .pipeline').not(showclasses.join(', ')).hide();
        if($('.pipelines-container .pipeline:visible').length == 0){ $('.no-pipelines').show(); }
        else { $('.no-pipelines').hide(); }
    });
    // Sort pipelines
    $('.pipelines-toolbar .pipeline-sorts button').click(function(){
        // Clicking sort a second time reverses the order
        var reverse = 1;
        if($(this).hasClass('active') && $(this).hasClass('reverse')){
            var reverse = -1;
            $(this).removeClass('reverse');
        } else {
            $('.pipelines-toolbar .pipeline-sorts button').removeClass('reverse');
            $(this).addClass('reverse');
        }
        // Apply the active class to this button
        $('.pipelines-toolbar .pipeline-sorts button').removeClass('active');
        $(this).blur().addClass('active');
        // Sort the pipeline cards
        $pipelines = $('.pipelines-container .pipeline');
        if($(this).text() == 'Alphabetical'){
            $pipelines.sort(function(a,b){
                var an = $(a).find('.card-title .pipeline-name').text();
                var bn = $(b).find('.card-title .pipeline-name').text();
                if(an > bn) { return 1 * reverse; }
                if(an < bn) { return -1 * reverse; }
                return 0;
            });
        }
        if($(this).text() == 'Status'){
            $pipelines.sort(function(a,b){
                var at = $(a).attr("class").match(/pipeline-[\w]*\b/)[0];
                var bt = $(b).attr("class").match(/pipeline-[\w]*\b/)[0];
                if(at == 'pipeline-released'){ an = 3; }
                if(at == 'pipeline-archived'){ an = 2; }
                if(at == 'pipeline-dev'){ an = 1; }
                if(bt == 'pipeline-released'){ bn = 3; }
                if(bt == 'pipeline-archived'){ bn = 2; }
                if(bt == 'pipeline-dev'){ bn = 1; }
                if(an > bn) { return -1 * reverse; }
                if(an < bn) { return 1 * reverse; }
                return 0;
            });
        }
        if($(this).text() == 'Stars'){
            $pipelines.sort(function(a,b){
                var an = Number($(a).find('.stargazers').text().trim());
                var bn = Number($(b).find('.stargazers').text().trim());
                if(typeof an === "undefined" || isNaN(an)){ an = 0; }
                if(typeof bn === "undefined" || isNaN(bn)){ bn = 0; }
                if(an > bn) { return -1 * reverse; }
                if(an < bn) { return 1 * reverse; }
                return 0;
            });
        }
        if($(this).text() == 'Last Release'){
            $pipelines.sort(function(a,b){
                var an = $(a).find('.publish-date').data('pubdate');
                var bn = $(b).find('.publish-date').data('pubdate');
                if(typeof an === "undefined"){ an = 0; }
                if(typeof bn === "undefined"){ bn = 0; }
                if(an > bn) { return -1 * reverse; }
                if(an < bn) { return 1 * reverse; }
                return 0;
            });
        }
        $pipelines.detach().appendTo($('.pipelines-container'));
    });
});
