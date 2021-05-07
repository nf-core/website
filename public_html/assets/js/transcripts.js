//
// transcripts.js
// Custom javascript for combining transcripts with embedded YouTube videos
//

// Global variables
var player;

    
      function onYouTubeIframeAPIReady() {
        player = new YT.Player("video-placeholder", {
          width: 600,
          height: 400,
          videoId: video_id,
          playerVars: {
            playsinline: 1,
            color: "green"
          },
          events: {
            onReady: initialize,
          },
        });
      }
    function initialize() {
    //    event.target.playVideo();
    
    // player.playVideo();
    $('details a[href^="https://youtu"]').each(function () {
      $(this).data("timestamp", $(this.href.split(";t=")).last()[0]);

      if ($(this).data("timestamp")) {
        $(this).addClass("timestamp-link");
        this.href = "javascript:void(0)";
      }
    });
    $("a.timestamp-link").click(function (e) {
      e.preventDefault();
      debugger;
      player.seekTo($(e).data("timestamp"), true);
    });
    }
