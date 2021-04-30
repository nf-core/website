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
    player.seekTo(200,true);
    player.playVideo();

    }
