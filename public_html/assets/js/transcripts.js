//
// transcripts.js
// Custom javascript for combining transcripts with embedded YouTube videos
//

// Global variables
var player;

    $('a[href^="https://youtu"]').each(function(a){
      debugger;
      a.data('timestamp',a.attr('href').split(";t=").last())
      
      if(a.data('timestamp')){
        a.addClass('timestamp-link')
      }
    });
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
    player.seekTo(200,true);
    player.playVideo();

    }
