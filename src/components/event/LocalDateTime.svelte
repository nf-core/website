<script lang="ts">
    export let date: Date;
    export let date_options: Intl.DateTimeFormatOptions = {
        year: 'numeric',
        month: 'long',
        day: 'numeric',
        hour: 'numeric',
        minute: 'numeric',
        hour12: false,
    };
    date = new Date(date + 'Z');
</script>

<span>{date.toLocaleString('en-US', date_options)}</span>
