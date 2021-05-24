//
// transcripts.js
// Custom javascript for combining transcripts with embedded YouTube videos
//

// Global variables
var player;
var timer;
var current_highlight;

    
function onYouTubeIframeAPIReady() {
  player = new YT.Player("video-placeholder", {
    width: 560,
    height: 315,
    videoId: video_id,
    playerVars: {
      playsinline: 1,
      color: "green",
    },
    events: {
      onReady: initialize,
      onStateChange: highlight,
    },
  });
}
function initialize() {
$('details a[href^="https://youtu"],details a[href^="https://www.youtu"]').each(function () {
  let timestamp = $(this.href.split(";t=")).last()[0];
  if (timestamp.indexOf("s") > 0 || timestamp.indexOf("m") > 0){
    timestamp =
      parseInt(timestamp.substring(0, timestamp.indexOf("m")) * 60) +
      parseInt(
        timestamp.substring(timestamp.indexOf("m") + 1, timestamp.indexOf("s"))
      );
  } else {
    timestamp = parseInt(timestamp)
  }
  $(this).data("timestamp", timestamp);

  if ($(this).data("timestamp")) {
    $(this).addClass("timestamp-link");
    this.href = "javascript:void(0)";
  }
});
$("a.timestamp-link").click(function (e) {
  e.preventDefault();
  player.seekTo($(e.target).data("timestamp"), true);
});
}
function highlight(event) {
  if (event.data === 1) {
    timer = setInterval(function () {
      var current_time = player.getCurrentTime();
      let all_highlights = $("a.timestamp-link").filter(function () {
        return $(this).data("timestamp") < current_time;
      });
      if ( current_highlight === undefined || all_highlights.last()[0]!== current_highlight[0]){
        $(current_highlight).parent("p").removeClass("bg-success-light");
        current_highlight = all_highlights.last();
        $(current_highlight).parent("p").addClass("bg-success-light");
      }
    }, 500);
  } else if(event.data===2) {
    clearInterval(timer);
  }
  
}