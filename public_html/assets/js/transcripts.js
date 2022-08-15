//
// transcripts.js
// Custom javascript for combining transcripts with embedded YouTube videos
//

// Global variables
var player;
var timer;
var current_highlight;

var tag = document.createElement('script');

tag.src = 'https://www.youtube.com/iframe_api';
var firstScriptTag = document.getElementsByTagName('script')[0];
firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);

function onYouTubeIframeAPIReady() {
  player = new YT.Player('video-placeholder', {
    width: 560,
    height: 315,
    videoId: video_id,
    playerVars: {
      playsinline: 1,
      color: 'green',
    },
    events: {
      onReady: initialize,
      onStateChange: highlight,
    },
  });
}
function initialize() {
  $('details a[href^="https://youtu"],details a[href^="https://www.youtu"]').each(function () {
    let timestamp = $(this.href.split(';t=')).last()[0];
    if (timestamp.indexOf('s') > 0 || timestamp.indexOf('m') > 0) {
      timestamp =
        parseInt(timestamp.substring(0, timestamp.indexOf('m')) * 60) +
        parseInt(timestamp.substring(timestamp.indexOf('m') + 1, timestamp.indexOf('s')));
    } else {
      timestamp = parseInt(timestamp);
    }
    $(this).data('timestamp', timestamp);

    if ($(this).data('timestamp')) {
      $(this).addClass('timestamp-link');
      this.href = 'javascript:void(0)';
    }
  });
  $('details').prop('open', true);
  $('a.timestamp-link').click(function (e) {
    e.preventDefault();
    player.seekTo($(e.target).data('timestamp'), true);
  });
}
function highlight(event) {
  if (event.data === 1) {
    timer = setInterval(function () {
      var current_time = player.getCurrentTime();
      let all_highlights = $('a.timestamp-link').filter(function () {
        return $(this).data('timestamp') < current_time;
      });
      if (current_highlight === undefined || all_highlights.last()[0] !== current_highlight[0]) {
        $('p.highlighted').parent('div').removeClass('bg-success-light');
        $('p.highlighted').removeClass('highlighted');

        current_highlight = all_highlights.last();
        $(current_highlight).parent('p').addClass('highlighted');
        $('details p.highlighted')
          .nextAll()
          .each(function (i, element) {
            if ($(element).has('a.timestamp-link').length > 0) {
              return false;
            }
            $(element).addClass('highlighted');
          })
          .end();
        $('.highlighted').wrapAll("<div class='bg-success-light' />");
        if ($('details[open]').length > 0) {
          $('details[open]').animate(
            {
              scrollTop:
                $('.bg-success-light').offset().top -
                $('details[open]').offset().top +
                $('details[open]').scrollTop() -
                1.2 * parseFloat(getComputedStyle(document.documentElement).fontSize), // convert rem into px
            },
            'slow'
          );
        }
      }
    }, 500);
  } else if (event.data === 2) {
    clearInterval(timer);
  }
}
