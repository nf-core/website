<script>
    import OngoingEvents from './OngoingEvents.svelte';
    import { getCollection } from 'astro:content';

    export let events = [];
    export let event_type_classes = {};
    export let event_type_icons = {};

    const currentEvents = events.filter((event) => {
        let time_window = 24 * 60 * 60 * 1000;
        const event_start_unix = event.data.start.getTime();
        const event_end_unix = event.data.end.getTime();

        // increase time window to a week for events longer than 5 hours
        if (event_end_unix - event_start_unix > 5 * 60 * 60 * 1000) {
            time_window = 7 * 24 * 60 * 60 * 1000;
        }
        if (event.data.start < new Date() && new Date() < event.data.end) {
            return false;
        }

        if (event_start_unix < new Date().getTime() + time_window && new Date().getTime() < event_end_unix) {
        // TODO: uncoment above line and remove below line, which is only for testing purposes
        // if (new Date().getTime() < event_end_unix ) {
            return true;
        }
    });
console.log('currentEvents',currentEvents);
</script>

{#each currentEvents as event}
    <div class="current-event text-center">
        <h5>{event.data.title}</h5>
        <p class="card-text">{event.data.subtitle}</p>
        <a href={event.url} class="btn btn-primary"> Read More </a>
    </div>
{/each}

<style lang="scss">
</style>
