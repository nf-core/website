<script lang="ts">
    export let date: Date;
    export let date_options: Intl.DateTimeFormatOptions = {
        year: "numeric",
        month: "long",
        day: "numeric",
        hour: "numeric",
        minute: "numeric",
        hour12: false,
    };
    date = new Date(date);
    // don't show time if it date includes "00:00:00"
    if (date.toISOString().includes("00:00:00")) {
        date_options.hour = undefined;
        date_options.minute = undefined;
    }
</script>

<span>{date.toLocaleString("en-US", date_options)}</span>
