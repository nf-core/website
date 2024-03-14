<script lang="ts">
    import { onMount } from 'svelte';
    import YouTubePlayer from 'youtube-player';

    export let id: string = '';
    let timer: any;
    let current_highlight: any;
    let noTranscript = true;
    onMount(async () => {
        // get details element

        const details = document.querySelector('details');
        if (details === null) {
            return;
        } else {
            noTranscript = false;
        }
        // add video placeholder before details
        const videoPlaceholder = document.createElement('div');
        videoPlaceholder.id = 'video-placeholder';
        details.parentNode.insertBefore(videoPlaceholder, details);

        // create player
        let player = YouTubePlayer('video-placeholder', {
            videoId: id,
        });

        const timestampLinks = Array.from(
            document.querySelectorAll('details a[href^="https://youtu"], details a[href^="https://www.youtu"]'),
        );
        timestampLinks.forEach((link) => {
            let timestamp = link.href.split('&t=')[1];
            if (timestamp.includes('s') || timestamp.includes('m')) {
                const minutes = parseInt(timestamp.substring(0, timestamp.indexOf('m')));
                const seconds = parseInt(timestamp.substring(timestamp.indexOf('m') + 1, timestamp.indexOf('s')));
                timestamp = minutes * 60 + seconds;
            } else {
                timestamp = parseInt(timestamp);
            }

            link.dataset.timestamp = timestamp;
            link.classList.add('timestamp-link');
            link.href = 'javascript:void(0)';
        });
        document.querySelector('details').open = true;

        const allTimestampLinks = Array.from(document.querySelectorAll('a.timestamp-link'));
        allTimestampLinks.forEach((link) => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                player.seekTo(parseInt(e.target.dataset.timestamp), true);
            });
        });
        player.on('stateChange', (event) => {
            if (event.data === 1) {
                timer = setInterval(async () => {
                    const current_time = await player.getCurrentTime();
                    const timestampLinks = Array.from(document.querySelectorAll('a.timestamp-link'));

                    const all_highlights = timestampLinks.filter((link) => {
                        return parseInt(link.dataset.timestamp) < current_time;
                    });

                    if (
                        current_highlight === undefined ||
                        all_highlights[all_highlights.length - 1] !== current_highlight
                    ) {
                        const highlightedElements = Array.from(document.querySelectorAll('p.highlighted'));
                        if (highlightedElements.length > 0) {
                            highlightedElements.forEach((highlighted) => {
                                highlighted.classList.remove('highlighted', 'bg-success-light', 'text-white');
                            });
                        }
                        current_highlight = all_highlights[all_highlights.length - 1];
                        current_highlight.parentElement.classList.add('highlighted', 'bg-success-light', 'text-white');
                    }
                }, 500);
            } else if (event.data === 2) {
                clearInterval(timer);
            }
        });
    });
</script>

{#if noTranscript}
    <slot />
{/if}
