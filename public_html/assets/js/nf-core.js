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
        $('.pipelines-toolbar .pipeline-sorts button').removeClass('active');
        $(this).blur().addClass('active');
        $pipelines = $('.pipelines-container .pipeline');
        if($(this).text() == 'Alphabetical'){
            $pipelines.sort(function(a,b){
            	var an = $(a).find('.card-title .pipeline-name').text();
            	var bn = $(b).find('.card-title .pipeline-name').text();
            	if(an > bn) { return 1; }
            	if(an < bn) { return -1; }
            	return 0;
            });
        }
        if($(this).text() == 'Status'){
            $pipelines.sort(function(a,b){
            	var an = $(a).attr("class").match(/pipeline-[\w]*\b/)[0];
            	var bn = $(b).attr("class").match(/pipeline-[\w]*\b/)[0];
                if(an == bn) { return 0; }
                if(an == 'pipeline-released' && bn != 'pipeline-released'){ return 1; }
                if(an == 'pipeline-dev' && bn != 'pipeline-archived'){ return 1; }
            	return -1;
            });
        }
        if($(this).text() == 'Last Release'){
            $pipelines.sort(function(a,b){
                var an = $(a).find('.publish-date').data('pubdate');
            	var bn = $(b).find('.publish-date').data('pubdate');
                if(typeof an === "undefined"){ an = 0; }
                if(typeof bn === "undefined"){ bn = 0; }
                if(an > bn) { return -1; }
                if(an < bn) { return 1; }
                return 0;
            });
        }
        $pipelines.detach().appendTo($('.pipelines-container'));
    });
})
