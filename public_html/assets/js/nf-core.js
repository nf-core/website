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

    // Function to switch CSS theme file
    $('.theme-switcher label').click(function(){
        var theme = $(this).find('input').val();
        var newlink = '/assets/css/nf-core-'+theme+'.css';
        $('#theme-stylesheet').attr('href', newlink);
        document.cookie = 'nfcoretheme='+theme+'; expires=Thu, 2 Dec 2032 12:00:00 UTC; path=/';
    });

    // Override the .contains() filter to be case insenstive
    $.expr[":"].contains = $.expr.createPseudo(function(arg) {
        return function( elem ) {
            return $(elem).text().toUpperCase().indexOf(arg.toUpperCase()) >= 0;
        };
    });

    // Homepage video switcher
    $('.video-chooser a').click(function(e){
        if($('#nf-core-video').is(':visible')){
          e.preventDefault();
          $('.video-chooser a').removeClass('active');
          $(this).addClass('active');
          $('#nf-core-video').attr('src', $(this).data('src'));
        }
    });

    // Homepage contributor images fading in and out
    var h_contrib_imgs = $('.homepage_contrib_logos a');
    if(h_contrib_imgs.length > 0){
        setTimeout(function(){ switch_contrib_img(); }, 2000);
    }
    function switch_contrib_img(){
        // Reset if all images have been shown
        if($('.homepage_contrib_logos a:hidden:not(.contrib_shown)').length == 0){
            $('.homepage_contrib_logos a').removeClass('contrib_shown');
            $('.homepage_contrib_logos a:visible').addClass('contrib_shown');
        }

        // Get random seeds for images to fade in and out
        var vis_imgs = $('.homepage_contrib_logos a:visible');
        var img_out_idx = Math.floor(Math.random()*vis_imgs.length);

        var hidden_imgs = $('.homepage_contrib_logos a:hidden:not(.contrib_shown)');
        var img_in_idx = Math.floor(Math.random()*hidden_imgs.length);

        // Label with a class
        hidden_imgs.eq(img_in_idx).addClass('contrib_fade_in contrib_shown');
        vis_imgs.eq(img_out_idx).addClass('contrib_fade_out');

        // Move image to be faded in next to the one to be faded out
        hidden_imgs.eq(img_in_idx).detach().insertAfter(vis_imgs.eq(img_out_idx));

        // Fade images in and out
        $('.contrib_fade_in').fadeIn();
        $('.contrib_fade_out').hide();

        // Clear labels
        $('.homepage_contrib_logos a').removeClass('contrib_fade_in contrib_fade_out');

        // Run this again in 2 seconds
        setTimeout(function(){ switch_contrib_img(); }, 3000);
    }

    // Filter pipelines with text
    function filter_pipelines_text(ftext){
        $('.pipelines-container .pipeline:contains("'+ftext+'")').show();
        $('.pipelines-container .pipeline:not(:contains("'+ftext+'"))').hide();
        if($('.pipelines-container .pipeline:visible').length == 0){ $('.no-pipelines').show(); }
        else { $('.no-pipelines').hide(); }
    }
    // page load
    if($('.pipelines-toolbar .pipeline-filters input').val()){
        var ftext = $('.pipelines-toolbar .pipeline-filters input').val();
        filter_pipelines_text(ftext);
        $('.pipelines-toolbar .pipeline-filters input').addClass('active');
    }
    // onchange
    $('.pipelines-toolbar .pipeline-filters input').keyup(function(){
        var ftext = $('.pipelines-toolbar .pipeline-filters input').val();
        filter_pipelines_text(ftext);
        if($('.pipelines-toolbar .pipeline-filters input').val()){
            $('.pipelines-toolbar .pipeline-filters input').addClass('active');
        } else {
            $('.pipelines-toolbar .pipeline-filters input').removeClass('active');
        }
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
    // Change pipelines display type
    $('.display-btn').click(function(){
        $('.display-btn').removeClass('active');
        $(this).addClass('active');
        var dtype = $(this).data('dtype');
        if(dtype == 'list' || dtype == 'blocks'){
            $('.pipelines-container').
                removeClass('pipelines-container-blocks pipelines-container-list').
                addClass('pipelines-container-'+dtype);
        }
    });

    // Make the stats tables sortable
    $('.pipeline-stats-table').tablesorter();
});
