//
// nf-core.js
// Custom javascript for the nf-core website
//

$(function () {
  // Enable tooltips
  $('.container').tooltip({
    //can't use body here, because scrollspy has already an event on it and bootstrap only allows one per selector.
    selector: '[data-bs-toggle="tooltip"]',
    html: true,
  });
  if (document.querySelector('.toc')) {
    // Enable scrollspy
    var scrollSpy = new bootstrap.ScrollSpy(document.body, {
      target: '.toc>div>div>.nav',
    });
  }
  // Enable code highlighting

  // Don't try to guess markdown language to highlight (gets it wrong most of the time)
  hljs.configure({ languages: [], ignoreUnescapedHTML: true });
  // don't highlight certain keywords
  const REMOVE_KEYWORDS = ['test', 'help'];
  hljs.getLanguage('bash').keywords.built_in = hljs.getLanguage('bash').keywords.built_in.filter((element) => {
    return !REMOVE_KEYWORDS.includes(element);
  });
  // do highlight custom keywords
  const CUSTOM_KEYWORDS = ['nextflow'];
  hljs.getLanguage('bash').keywords.built_in.push(...CUSTOM_KEYWORDS);

  hljs.highlightAll();

  // add copy-paste button to code blocks
  $('.hljs.language-bash').each(function () {
    let $this = $(this);
    $this.addClass('position-relative');
    var $copy_btn = $(
      '<button class="btn btn-outline-secondary copy-button position-absolute top-0 end-0" data-bs-toggle="tooltip"  data-bs-placement="left" data-bs-original-title="Copy to clipboard" ><i class="fas fa-clipboard fa-swap-opacity"></i></button>',
    );
    $this.append($copy_btn);
    var btn_tooltip = new bootstrap.Tooltip($copy_btn, { container: 'body' });
    $copy_btn.on('click', function () {
      let $text = $(this).parent('code.language-bash').text();
      navigator.clipboard.writeText($text).then(() => {
        // debugger;
        btn_tooltip._element.setAttribute('data-bs-original-title', 'Copied!');
        btn_tooltip.show();
        // reset the tooltip title
        btn_tooltip._element.setAttribute('data-bs-original-title', 'Copy to clipboard');
      });
    });
  });

  // Set theme cookie if not set
  if (document.cookie.indexOf('nfcoretheme') == -1) {
    if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
      document.cookie = 'nfcoretheme=dark; expires=Thu, 2 Dec 2032 12:00:00 UTC; path=/';
      // Run the function so that dark-mode images are switched in
      update_theme('dark');
    } else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) {
      document.cookie = 'nfcoretheme=light; expires=Thu, 2 Dec 2032 12:00:00 UTC; path=/';
    }
  }
  // update cookie when OS theme changes
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', (e) => {
    var new_theme = e.matches ? 'dark' : 'light';
    update_theme(new_theme);
  });

  $('.theme-switcher label').on('click', function () {
    var theme = $(this).attr('for').split('-')[1];

    //uncheck all radio buttons and select only current one
    $('.theme-switcher input:checked').prop('checked', false);
    $('.theme-switcher #theme-' + theme).prop('checked', true);
    update_theme(theme);
  });
  // Function to switch CSS theme
  function update_theme(theme) {
    // Switch the stylesheet
    var newlink = '/assets/css/nf-core-' + theme + '.css';
    $('#theme-stylesheet').attr('href', newlink);

    // Switch any images
    $('.darkmode-image').each(function () {
      var old_imgdir = theme == 'dark' ? 'contributors-colour' : 'contributors-white';
      var new_imgdir = theme == 'dark' ? 'contributors-white' : 'contributors-colour';
      $(this).attr('src', $(this).attr('src').replace(old_imgdir, new_imgdir));
    });

    // Set a cookie to remember
    document.cookie = 'nfcoretheme=' + theme + '; expires=Thu, 2 Dec 2032 12:00:00 UTC; path=/';
  }
  // Override the .contains() filter to be case insensitive
  $.expr[':'].contains = $.expr.createPseudo(function (arg) {
    return function (elem) {
      return $(elem).text().toUpperCase().indexOf(arg.toUpperCase()) >= 0;
    };
  });

  // Homepage video switcher
  $('.video-chooser a').click(function (e) {
    if ($('#nf-core-video').is(':visible')) {
      e.preventDefault();
      $('.video-chooser a').removeClass('active');
      $(this).addClass('active');
      $('#nf-core-video').attr('src', $(this).data('src'));
    }
  });

  // Homepage contributor images fading in and out
  var h_contrib_imgs = $('.homepage_contrib_logos a');
  if (h_contrib_imgs.length > 0) {
    setTimeout(function () {
      switch_contrib_img();
    }, 2000);
  }
  function switch_contrib_img() {
    // Add image that will be removed to the end of the list
    contributors_imgs.push($('.homepage_contrib_logos a:first-child').html());

    // Animate the first image off the screen and then remove
    var margin_left = $('.homepage_contrib_logos a:first-child').width() * -1;
    $('.homepage_contrib_logos a:first-child').animate({ 'margin-left': margin_left }, function () {
      $(this).remove();
    });

    // Add a new image to the end of the list (should be off the screen to the right)
    var next_img = contributors_imgs.shift();
    $('.homepage_contrib_logos').append($(next_img));

    // Run this again in 2 seconds
    setTimeout(function () {
      switch_contrib_img();
    }, 3000);
  }

  // Filter pipelines with text
  function filter_pipelines_text(ftext) {
    $('.pipelines-container .pipeline:contains("' + ftext + '")')
      .parent('.col')
      .show();
    $('.pipelines-container .pipeline:not(:contains("' + ftext + '"))')
      .parent('.col')
      .hide();
    if ($('.pipelines-container .pipeline:visible').length == 0) {
      $('.no-pipelines').show();
    } else {
      $('.no-pipelines').hide();
    }
  }
  // page load
  if ($('.pipelines-toolbar .pipeline-filters input').val()) {
    var ftext = $('.pipelines-toolbar .pipeline-filters input').val().trim();
    filter_pipelines_text(ftext);
    $('.pipelines-toolbar .pipeline-filters input').addClass('active');
  }
  // onchange
  $('.pipelines-toolbar .pipeline-filters input').keyup(function () {
    var ftext = $('.pipelines-toolbar .pipeline-filters input').val().trim();
    filter_pipelines_text(ftext);
    if ($('.pipelines-toolbar .pipeline-filters input').val()) {
      $('.pipelines-toolbar .pipeline-filters input').addClass('active');
    } else {
      $('.pipelines-toolbar .pipeline-filters input').removeClass('active');
    }
  });

  // Filter pipelines with buttons
  $('.pipelines-toolbar .pipeline-filters button').click(function () {
    $(this).blur().toggleClass('active');
    var showclasses = [];
    $('.pipelines-toolbar .pipeline-filters button.active').each(function () {
      showclasses.push($(this).data('bsTarget'));
    });
    $('.pipelines-container .pipeline').filter(showclasses.join(', ')).parent('.col').show();
    $('.pipelines-container .pipeline').not(showclasses.join(', ')).parent('.col').hide();
    if ($('.pipelines-container .pipeline:visible').length == 0) {
      $('.no-pipelines').show();
    } else {
      $('.no-pipelines').hide();
    }
  });
  // Sort pipelines
  $('.pipelines-toolbar .pipeline-sorts button').click(function () {
    // Clicking sort a second time reverses the order
    var reverse = 1;
    if ($(this).hasClass('active') && $(this).hasClass('reverse')) {
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
    $pipelines = $('.pipelines-container .col');
    if ($(this).text() == 'Alphabetical') {
      $pipelines.sort(function (a, b) {
        var an = $(a).find('.card-title .pipeline-name').text();
        var bn = $(b).find('.card-title .pipeline-name').text();
        if (an > bn) {
          return 1 * reverse;
        }
        if (an < bn) {
          return -1 * reverse;
        }
        return 0;
      });
    }
    if ($(this).text() == 'Status') {
      $pipelines.sort(function (a, b) {
        var at = $(a)
          .attr('class')
          .match(/pipeline-[\w]*\b/)[0];
        var bt = $(b)
          .attr('class')
          .match(/pipeline-[\w]*\b/)[0];
        if (at == 'pipeline-released') {
          an = 3;
        }
        if (at == 'pipeline-archived') {
          an = 2;
        }
        if (at == 'pipeline-dev') {
          an = 1;
        }
        if (bt == 'pipeline-released') {
          bn = 3;
        }
        if (bt == 'pipeline-archived') {
          bn = 2;
        }
        if (bt == 'pipeline-dev') {
          bn = 1;
        }
        if (an > bn) {
          return -1 * reverse;
        }
        if (an < bn) {
          return 1 * reverse;
        }
        return 0;
      });
    }
    if ($(this).text() == 'Stars') {
      $pipelines.sort(function (a, b) {
        var an = Number($(a).find('.stargazers').text().trim());
        var bn = Number($(b).find('.stargazers').text().trim());
        if (typeof an === 'undefined' || isNaN(an)) {
          an = 0;
        }
        if (typeof bn === 'undefined' || isNaN(bn)) {
          bn = 0;
        }
        if (an > bn) {
          return -1 * reverse;
        }
        if (an < bn) {
          return 1 * reverse;
        }
        return 0;
      });
    }
    if ($(this).text() == 'Last Release') {
      $pipelines.sort(function (a, b) {
        var an = $(a).find('.publish-date').data('pubdate');
        var bn = $(b).find('.publish-date').data('pubdate');
        if (typeof an === 'undefined') {
          an = 0;
        }
        if (typeof bn === 'undefined') {
          bn = 0;
        }
        if (an > bn) {
          return -1 * reverse;
        }
        if (an < bn) {
          return 1 * reverse;
        }
        return 0;
      });
    }
    $pipelines.detach().appendTo($('.pipelines-container'));
  });
  // Change pipelines display type
  $('.display-btn').click(function () {
    $('.display-btn').removeClass('active');
    $(this).addClass('active');
    var dtype = $(this).data('dtype');
    if (dtype == 'list' || dtype == 'blocks') {
      $('.pipelines-container')
        .removeClass('pipelines-container-blocks pipelines-container-list')
        .addClass('pipelines-container-' + dtype)
        .removeClass('row-cols-md-2 g-4');
      if (dtype === 'blocks') {
        $('.pipelines-container').addClass('row-cols-md-2 g-4');
      }
    }
  });
  // Filter modules with text
  function filter_modules_text(ftext) {
    $('.modules-container .module:contains("' + ftext + '")').show();
    $('.modules-container .module:not(:contains("' + ftext + '"))').hide();
    if ($('.modules-container .module:visible').length == 0) {
      $('.no-modules').show();
    } else {
      $('.no-modules').hide();
    }
    update_facet_bar();
  }
  // Filter modules with facets
  function filter_modules_facet(ftext) {
    // undo filtering if repeated click or no text
    if (ftext.length == 0) {
      $('.modules-container .card').show();
    } else {
      $('.modules-container .module .topics').each(function () {
        $('.badge', this).each(function () {
          if (ftext.includes($(this).text())) {
            $(this).parents('.card').show();
            return false;
          } else $(this).parents('.card').hide();
        });
      });
    }
    update_facet_bar();
    return true;
  }

  function update_facet_bar() {
    var facets = [];
    $('.modules-container .module:visible .topics .badge').each(function () {
      var facet = $(this).text();
      facets.push(facet);
    });
    var facet_occurrences = facets.reduce((fac, cur) => {
      fac[cur] ??= 0;
      fac[cur]++;
      return fac;
    }, {});
    $('.facet-bar .facet-item').each(function () {
      if (Object.keys(facet_occurrences).includes($('.facet-name', this).text())) {
        $('.facet-value', this).text(facet_occurrences[$('.facet-name', this).text()]);
      } else {
        $('.facet-value', this).text('0');
      }
    });
    //sort facets by count
    $('.facet-bar .facet-item')
      .sort(function (a, b) {
        return $('.facet-value', b).text() - $('.facet-value', a).text();
      })
      .appendTo('.facet-bar ul');
    return true;
  }
  // page load
  if ($('.modules-toolbar .module-filters input').val()) {
    var ftext = $('.modules-toolbar .module-filters input').val().trim();
    filter_modules_text(ftext);
    $('.modules-toolbar .module-filters input').addClass('active');
  }
  // onchange
  $('.modules-toolbar .module-filters input').on('keyup', function () {
    var ftext = $('.modules-toolbar .module-filters input').val().trim();
    filter_modules_text(ftext);
  });
  // on facet click
  $('.facet-bar .facet-item').on('click', function (e) {
    // undo selection on click on active items
    if ($(e.currentTarget).hasClass('active')) {
      $(e.currentTarget).removeClass('active');
    } else {
      $(e.currentTarget).addClass('active');
    }
    // join text from active facets with semi-colon
    var ftext = $('.facet-bar .facet-item.active .facet-name')
      .map(function () {
        return $(this).text();
      })
      .get();

    filter_modules_facet(ftext);
  });

  // on badge click
  $('.topics .badge').on('click', function (e) {
    var keyword = $(e.currentTarget).data('keyword');
    $('.facet-bar #' + keyword).trigger('click');
  });

  // Make the stats tables sortable
  if ($('.pipeline-stats-table').length > 0) {
    $('.pipeline-stats-table').tablesorter();
  }

  // Pipeline page version number dropdown
  $('#version_select').on('change', function () {
    document.location.href = $(this).val();
  });

  //copy text to clipboard
  $('.toast').toast();
  $('.copy-txt').on('click', function () {
    var target = $(this).data('bsTarget');
    var target_id = '#' + target;
    $(target_id).trigger('select');
    document.execCommand('copy');
    $(target_id).trigger('blur');
    $('.cmd_copied').toast('show');
  });
  if (window.location.hash & ($('.param-docs').length > 0)) {
    scroll_to($(window.location.hash), 0);
  }
  // Page-scroll links
  $('body').on('click', '.scroll_to_link', function (e) {
    e.preventDefault();
    var current_href = $(this).attr('href');
    history.replaceState(null, '', current_href); // add href to url
    scroll_to($(current_href), 0);
  });

  if ($('details>summary:contains("Video transcript")').length > 0) {
    $('details>summary:contains("Video transcript")')
      .parents('details')
      .before('<div class="ratio ratio-16x9"><div id="video-placeholder"></div></div>');
  } else if ($('.rendered-markdown').length > 0 && typeof youtube_embed !== 'undefined' && youtube_embed) {
    $('.rendered-markdown').append('<div class="ratio ratio-16x9"><div id="video-placeholder"></div></div>');
  }

  // Expand all details on page
  $('.expand-details').on('click', function () {
    // if all details are already open, close them, else open all
    if ($('details[open]').length === $('details').length) {
      $('details[open]').removeAttr('open');
    } else {
      $('details:not([open])').attr('open', '');
    }
    // refresh scrollspy to recalculate offsets, still not working 100% of the time...
    var dataSpyList = [].slice.call(document.querySelectorAll('[data-bs-spy="scroll"]'));
    dataSpyList.forEach(function (dataSpyEl) {
      bootstrap.ScrollSpy.getInstance(dataSpyEl).refresh();
    });
  });

  // Filter events with buttons
  $('.events-toolbar .events-filters input[type="radio"]').click(function () {
    $(this).blur().toggleClass('active');
    var showclasses = $(this).data('bsTarget') ? $(this).data('bsTarget') : [];
    if (showclasses.length > 0) {
      $('.event-list .card').filter(showclasses).show();
      $('.event-list .card').not(showclasses).hide();
    } else {
      // reset filtering
      $('.event-list .card').show();
    }
  });

  // Filter pipelines with text
  function filter_events_text(ftext) {
    if (ftext.length > 0) {
      $('.event-list .card:visible:contains("' + ftext + '")').show();
      $('.event-list .card:visible:not(:contains("' + ftext + '"))').hide();
    } else {
      // reset filtering
      $('.event-list .card').show();
    }
  }
  // page load
  if ($('.events-toolbar .event-filters input').val()) {
    var ftext = $('.events-toolbar .event-filters input').val();
    filter_events_text(ftext);
    $('.events-toolbar .event-filters input').addClass('active');
  }
  // onchange
  $('.events-toolbar .event-filters input').keyup(function () {
    var ftext = $('.events-toolbar .event-filters input').val();
    filter_events_text(ftext);
    if ($('.events-toolbar .event-filters input').val()) {
      $('.events-toolbar .event-filters input').addClass('active');
    } else {
      $('.events-toolbar .event-filters input').removeClass('active');
    }
  });
});

function scroll_to(target_el, offset) {
  var el_offset = target_el.offset().top - offset;
  $([document.documentElement, document.body]).animate(
    {
      scrollTop: el_offset,
    },
    500,
  );
}
