<script lang="ts">
    interface Props {
        date: Date;
        date_options?: Intl.DateTimeFormatOptions;
        showTimezone?: boolean;
    }

    let {
        date = $bindable(),
        date_options = {
            year: "numeric",
            month: "long",
            day: "numeric",
            hour: "numeric",
            minute: "numeric",
            hour12: false,
        },
        showTimezone = false,
    }: Props = $props();

    // Get user's local timezone
    const userTimezone = Intl.DateTimeFormat().resolvedOptions().timeZone;

    // don't show time and don't convert to local time if it date includes "00:00:00"
    const formattedOptions = $derived(
        date.toISOString().includes("00:00:00")
            ? { ...date_options, hour: undefined, minute: undefined, timeZone: "UTC" }
            : date_options,
    );
</script>

<span>
    {date.toLocaleString("en-US", formattedOptions)}
    {#if showTimezone}
        in {userTimezone}
    {/if}
</span>
